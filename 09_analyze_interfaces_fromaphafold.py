"""
For each species (relaxed_model_*.pdb):
  - Buried surface area for all 6 protein pairs (AB, AC, AD, BC, BD, CD)
  - Interface residues for each pair
  - Contact classification: H-bonds, salt bridges, hydrophobic, van der Waals

Outputs per species (written to <species_dir>/analysis_py/):
  - interface_areas.txt
  - contacts_{pair}.tsv      (atom-level contacts)
  - residues_{pair}.tsv      (residue-level summary)
  - summary_{species}.txt    (human-readable summary)

Outputs overall (written to this script's directory):
  - summary_interface_areas.csv
  - summary_contact_counts.csv
"""

import os, sys, copy, warnings, glob, csv, math
from collections import defaultdict

warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# ─── constants ────────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

PROTEIN_PAIRS = [("A","B"),("A","C"),("A","D"),("B","C"),("B","D"),("C","D")]

HBOND_POLAR = {"N", "O", "S"}

POS_CHARGED = {
    "ARG": {"NH1","NH2","NE","CZ"},
    "LYS": {"NZ"},
    "HIS": {"ND1","NE2"},
}
NEG_CHARGED = {
    "ASP": {"OD1","OD2"},
    "GLU": {"OE1","OE2"},
}

HYDROPHOBIC_RES = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","TYR","GLY"}

CONTACT_CUTOFF     = 5.0   # Å – general vdW/contact envelope
CLASH_CUTOFF       = 2.0   # Å – below this = steric clash
HBOND_CUTOFF       = 3.5   # Å – heavy-atom H-bond donor…acceptor
SALTBRIDGE_CUTOFF  = 5.5   # Å – charged atom distance
HYDROPHOBIC_CUTOFF = 4.5   # Å – C…C hydrophobic


# ─── helpers ──────────────────────────────────────────────────────────────────

def atom_dist(a1, a2):
    diff = a1.coord - a2.coord
    return math.sqrt(float(diff.dot(diff)))


def build_chain_atoms(model, chain_ids):
    """Return list of (chain_id, residue, atom) for given chains, no hydrogens."""
    atoms = []
    for ch in model.get_chains():
        if ch.id in chain_ids:
            for res in ch.get_residues():
                for atom in res.get_atoms():
                    if atom.element == "H" or atom.name.startswith("H"):
                        continue
                    atoms.append((ch.id, res, atom))
    return atoms


def copy_sub(model, chain_list):
    mc = copy.deepcopy(model)
    for c in list(mc.get_chains()):
        if c.id not in chain_list:
            mc.detach_child(c.id)
    return mc


def calc_buried_area(full_model, chains1, chains2):
    """Buried surface area = (SASA_g1 + SASA_g2 - SASA_g1g2) / 2."""
    sr = ShrakeRupley()

    m12 = copy_sub(full_model, chains1 + chains2)
    m1  = copy_sub(full_model, chains1)
    m2  = copy_sub(full_model, chains2)

    sr.compute(m12, level="A")
    sasa_12 = sum(a.sasa for c in m12 for r in c for a in r)

    sr.compute(m1, level="A")
    sasa_1 = sum(a.sasa for c in m1 for r in c for a in r)

    sr.compute(m2, level="A")
    sasa_2 = sum(a.sasa for c in m2 for r in c for a in r)

    return (sasa_1 + sasa_2 - sasa_12) / 2.0


def get_interface_residues(full_model, chains1, chains2, delta_sasa_cutoff=1.0):
    """Return residues whose SASA decreases by >= cutoff upon complex formation."""
    sr = ShrakeRupley()

    m12 = copy_sub(full_model, chains1 + chains2)
    m1  = copy_sub(full_model, chains1)
    m2  = copy_sub(full_model, chains2)

    sr.compute(m12, level="R")
    sasa_complex = {(c.id, r.id[1]): r.sasa for c in m12 for r in c}

    sr.compute(m1, level="R")
    sasa_g1 = {(c.id, r.id[1]): r.sasa for c in m1 for r in c}

    sr.compute(m2, level="R")
    sasa_g2 = {(c.id, r.id[1]): r.sasa for c in m2 for r in c}

    interface1, interface2 = [], []

    for c in m1.get_chains():
        for r in c.get_residues():
            key = (c.id, r.id[1])
            delta = sasa_g1.get(key, 0) - sasa_complex.get(key, 0)
            if delta >= delta_sasa_cutoff:
                interface1.append((c.id, r.get_resname(), r.id[1], round(delta, 2)))

    for c in m2.get_chains():
        for r in c.get_residues():
            key = (c.id, r.id[1])
            delta = sasa_g2.get(key, 0) - sasa_complex.get(key, 0)
            if delta >= delta_sasa_cutoff:
                interface2.append((c.id, r.get_resname(), r.id[1], round(delta, 2)))

    return interface1, interface2


def classify_contact(res1, atom1, res2, atom2, dist):
    """
    Return interaction type string for a pair of atoms.
    Priority: clash > H-bond > salt bridge > hydrophobic > vdW
    """
    name1, name2 = atom1.name, atom2.name
    resn1, resn2 = res1.get_resname(), res2.get_resname()
    el1 = atom1.element if atom1.element else name1[0]
    el2 = atom2.element if atom2.element else name2[0]

    if dist < CLASH_CUTOFF:
        return "clash"

    if dist <= HBOND_CUTOFF and el1 in HBOND_POLAR and el2 in HBOND_POLAR:
        return "H-bond"

    is_pos1 = resn1 in POS_CHARGED and name1 in POS_CHARGED.get(resn1, set())
    is_neg1 = resn1 in NEG_CHARGED and name1 in NEG_CHARGED.get(resn1, set())
    is_pos2 = resn2 in POS_CHARGED and name2 in POS_CHARGED.get(resn2, set())
    is_neg2 = resn2 in NEG_CHARGED and name2 in NEG_CHARGED.get(resn2, set())

    if dist <= SALTBRIDGE_CUTOFF:
        if (is_pos1 and is_neg2) or (is_pos2 and is_neg1):
            return "salt bridge"

    if (el1 == "C" and el2 == "C" and
            resn1 in HYDROPHOBIC_RES and resn2 in HYDROPHOBIC_RES and
            dist <= HYDROPHOBIC_CUTOFF):
        return "hydrophobic"

    return "vdW"


def find_contacts(model, chains1, chains2):
    """Find all inter-group contacts within CONTACT_CUTOFF."""
    atoms1 = build_chain_atoms(model, chains1)
    atoms2 = build_chain_atoms(model, chains2)

    contacts = []
    for ch1, res1, atom1 in atoms1:
        for ch2, res2, atom2 in atoms2:
            d = atom_dist(atom1, atom2)
            if d > CONTACT_CUTOFF:
                continue
            ctype = classify_contact(res1, atom1, res2, atom2, d)
            contacts.append({
                "chain1":   ch1,
                "resname1": res1.get_resname(),
                "resnum1":  res1.id[1],
                "atom1":    atom1.name,
                "chain2":   ch2,
                "resname2": res2.get_resname(),
                "resnum2":  res2.id[1],
                "atom2":    atom2.name,
                "distance": round(d, 3),
                "type":     ctype,
            })
    return contacts


def deduplicate_contacts(contacts):
    """Remove duplicate atom pairs (keep shortest distance)."""
    seen = {}
    for c in contacts:
        key = tuple(sorted([
            (c["chain1"], c["resnum1"], c["atom1"]),
            (c["chain2"], c["resnum2"], c["atom2"])
        ]))
        if key not in seen or c["distance"] < seen[key]["distance"]:
            seen[key] = c
    return list(seen.values())


def contacts_to_residue_level(contacts):
    res_pairs = defaultdict(lambda: defaultdict(int))
    for c in contacts:
        key = (c["chain1"], c["resname1"], c["resnum1"],
               c["chain2"], c["resname2"], c["resnum2"])
        res_pairs[key][c["type"]] += 1
    return res_pairs


def count_by_type(contacts):
    counts = defaultdict(int)
    for c in contacts:
        counts[c["type"]] += 1
    return dict(counts)


# ─── writers ──────────────────────────────────────────────────────────────────

def write_contacts_tsv(contacts, filepath):
    if not contacts:
        with open(filepath, "w") as f:
            f.write("No contacts found\n")
        return
    fieldnames = ["chain1","resname1","resnum1","atom1",
                  "chain2","resname2","resnum2","atom2",
                  "distance","type"]
    type_order = {"clash":-1,"H-bond":0,"salt bridge":1,"hydrophobic":2,"vdW":3}
    contacts_sorted = sorted(contacts,
                             key=lambda c: (type_order.get(c["type"], 9), c["distance"]))
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(contacts_sorted)


def write_residue_contacts(contacts, iface1, iface2, filepath):
    res_pairs = contacts_to_residue_level(contacts)
    with open(filepath, "w") as f:
        f.write("=== Interface residues (group 1) ===\n")
        f.write(f"{'Chain':<6} {'Residue':<6} {'ResNum':<8} {'ΔSASA(Å²)':<12}\n")
        for ch, resn, resi, dsasa in sorted(iface1, key=lambda x: (x[0], x[2])):
            f.write(f"{ch:<6} {resn:<6} {resi:<8} {dsasa:<12}\n")

        f.write("\n=== Interface residues (group 2) ===\n")
        f.write(f"{'Chain':<6} {'Residue':<6} {'ResNum':<8} {'ΔSASA(Å²)':<12}\n")
        for ch, resn, resi, dsasa in sorted(iface2, key=lambda x: (x[0], x[2])):
            f.write(f"{ch:<6} {resn:<6} {resi:<8} {dsasa:<12}\n")

        f.write("\n=== Residue-pair contacts ===\n")
        header = (f"{'Chain1':<7}{'Res1':<6}{'Num1':<7}"
                  f"{'Chain2':<7}{'Res2':<6}{'Num2':<7}"
                  f"{'H-bond':<8}{'Salt-br':<9}{'Hydroph':<9}{'vdW':<6}")
        f.write(header + "\n")
        f.write("-"*65 + "\n")
        for (c1,rn1,ri1,c2,rn2,ri2), type_counts in sorted(
                res_pairs.items(), key=lambda x: (x[0][0], x[0][2], x[0][3], x[0][5])):
            line = (f"{c1:<7}{rn1:<6}{ri1:<7}"
                    f"{c2:<7}{rn2:<6}{ri2:<7}"
                    f"{type_counts.get('H-bond',0):<8}"
                    f"{type_counts.get('salt bridge',0):<9}"
                    f"{type_counts.get('hydrophobic',0):<9}"
                    f"{type_counts.get('vdW',0):<6}")
            f.write(line + "\n")


def write_species_summary(species_name, results, contact_summary, out_dir):
    outfile = os.path.join(out_dir, f"summary_{species_name}.txt")
    with open(outfile, "w") as f:
        f.write(f"Interface analysis summary: {species_name}\n")
        f.write("="*60 + "\n\n")

        f.write("BURIED SURFACE AREAS\n")
        f.write(f"{'Interface pair':<18} {'Buried area (Å²)':>18}\n")
        f.write("-"*38 + "\n")
        for pair, area in results["buried_areas"].items():
            f.write(f"{pair:<18} {area:>18.1f}\n")

        f.write("\n\nCONTACT COUNTS BY TYPE\n")
        f.write(f"{'Pair':<10} {'Clash':>7} {'H-bond':>8} "
                f"{'Salt-br':>8} {'Hydroph':>8} {'vdW':>8} {'Total':>8}\n")
        f.write("-"*62 + "\n")
        for label, contacts, iface1, iface2, cnt in contact_summary:
            total = sum(cnt.values())
            f.write(f"{label:<10} "
                    f"{cnt.get('clash',0):>7} "
                    f"{cnt.get('H-bond',0):>8} "
                    f"{cnt.get('salt bridge',0):>8} "
                    f"{cnt.get('hydrophobic',0):>8} "
                    f"{cnt.get('vdW',0):>8} "
                    f"{total:>8}\n")

        f.write("\n\nINTERFACE RESIDUES (delta-SASA >= 1 Å²)\n")
        for label, contacts, iface1, iface2, cnt in contact_summary:
            f.write(f"\n--- {label} ---\n")
            f.write("Group 1: " + ", ".join(
                f"{ch}/{resn}{resi}" for ch, resn, resi, _ in
                sorted(iface1, key=lambda x: (x[0], x[2]))
            ) + "\n")
            f.write("Group 2: " + ", ".join(
                f"{ch}/{resn}{resi}" for ch, resn, resi, _ in
                sorted(iface2, key=lambda x: (x[0], x[2]))
            ) + "\n")


# ─── per-species analysis ─────────────────────────────────────────────────────

def analyze_species(species_dir, species_name):
    """Run full interface analysis for one species."""
    pdb_files = sorted(glob.glob(os.path.join(species_dir, "relaxed_model_*.pdb")))
    if not pdb_files:
        print(f"  [SKIP] No relaxed_model_*.pdb found in {species_dir}")
        return None

    pdb_path = pdb_files[0]
    print(f"  Parsing {os.path.basename(pdb_path)} ...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(species_name, pdb_path)
    model = structure[0]

    available_chains = set(c.id for c in model.get_chains())
    print(f"  Chains: {sorted(available_chains)}")

    results = {"species": species_name, "buried_areas": {}, "contacts": {}}

    out_dir = os.path.join(species_dir, "analysis_py")
    os.makedirs(out_dir, exist_ok=True)

    # ── buried areas ──────────────────────────────────────────────────────────
    print("  Computing buried areas ...")
    area_lines = [f"Interface buried area analysis: {species_name}\n",
                  f"{'Pair':<10} {'Buried area (Å²)':>18}\n",
                  "-"*30 + "\n"]

    for ch1, ch2 in PROTEIN_PAIRS:
        if ch1 in available_chains and ch2 in available_chains:
            ba = calc_buried_area(model, [ch1], [ch2])
            label = f"{ch1}-{ch2}"
            results["buried_areas"][label] = round(ba, 1)
            area_lines.append(f"{label:<10} {ba:>18.1f}\n")
            print(f"    {label}: {ba:.1f} Å²")

    with open(os.path.join(out_dir, "interface_areas.txt"), "w") as f:
        f.writelines(area_lines)

    # ── contacts & interface residues ─────────────────────────────────────────
    print("  Computing protein-protein contacts ...")
    contact_summary = []

    for ch1, ch2 in PROTEIN_PAIRS:
        if ch1 not in available_chains or ch2 not in available_chains:
            continue
        label = f"{ch1}-{ch2}"
        print(f"    {label} ...")

        contacts = find_contacts(model, [ch1], [ch2])
        contacts = deduplicate_contacts(contacts)

        iface1, iface2 = get_interface_residues(model, [ch1], [ch2])

        write_contacts_tsv(contacts, os.path.join(out_dir, f"contacts_{label}.tsv"))
        write_residue_contacts(contacts, iface1, iface2,
                               os.path.join(out_dir, f"residues_{label}.tsv"))

        cnt = count_by_type(contacts)
        results["contacts"][label] = cnt
        contact_summary.append((label, contacts, iface1, iface2, cnt))

    write_species_summary(species_name, results, contact_summary, out_dir)
    return results


# ─── master summary ───────────────────────────────────────────────────────────

def write_master_csv(all_results, base_dir):
    area_rows = [r for r in all_results if r is not None]
    if area_rows:
        all_pairs = []
        for r in area_rows:
            for k in r["buried_areas"]:
                if k not in all_pairs:
                    all_pairs.append(k)
        fieldnames = ["species"] + all_pairs
        master_area = os.path.join(base_dir, "summary_interface_areas.csv")
        with open(master_area, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for r in area_rows:
                row = {"species": r["species"]}
                row.update(r["buried_areas"])
                writer.writerow(row)
        print(f"\nWrote {master_area}")

    contact_rows = []
    for r in all_results:
        if r is None:
            continue
        for pair, cnt in r["contacts"].items():
            contact_rows.append({
                "species":    r["species"],
                "pair":       pair,
                "clash":      cnt.get("clash", 0),
                "H-bond":     cnt.get("H-bond", 0),
                "salt_bridge":cnt.get("salt bridge", 0),
                "hydrophobic":cnt.get("hydrophobic", 0),
                "vdW":        cnt.get("vdW", 0),
                "total":      sum(cnt.values()),
            })

    if contact_rows:
        master_contacts = os.path.join(base_dir, "summary_contact_counts.csv")
        with open(master_contacts, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["species","pair","clash","H-bond",
                               "salt_bridge","hydrophobic","vdW","total"])
            writer.writeheader()
            writer.writerows(contact_rows)
        print(f"Wrote {master_contacts}")


# ─── main ─────────────────────────────────────────────────────────────────────

SPECIES_DIRS = [
    "arabidopsis", "artocarpus", "eleagnus_1", "humulus_1", "malus_1",
    "morus_1", "rosa", "ulmus_1", "urtica", "vitis", "zizyphus",
]

def main():
    species_dirs = []
    for name in SPECIES_DIRS:
        d = os.path.join(BASE_DIR, name)
        if os.path.isdir(d):
            species_dirs.append((d, name))

    # Also auto-discover any directories not in the explicit list
    for entry in sorted(os.listdir(BASE_DIR)):
        full = os.path.join(BASE_DIR, entry)
        if os.path.isdir(full) and entry not in SPECIES_DIRS and entry not in ("figures", "__pycache__", "snapshots", "chimeraX_results_hbond", "chimeraXcommands"):
            if glob.glob(os.path.join(full, "relaxed_model_*.pdb")):
                species_dirs.append((full, entry))

    if not species_dirs:
        print("No species directories with PDB files found.")
        sys.exit(1)

    print(f"Found {len(species_dirs)} species directories.\n")

    all_results = []
    for sp_dir, sp_name in species_dirs:
        print(f"\n{'='*60}")
        print(f"Species: {sp_name}")
        print(f"{'='*60}")
        result = analyze_species(sp_dir, sp_name)
        all_results.append(result)

    print(f"\n{'='*60}")
    print("Writing master summary files ...")
    write_master_csv(all_results, BASE_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
