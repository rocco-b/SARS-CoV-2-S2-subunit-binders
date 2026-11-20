import os

def split_pdb(file_path, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(file_path, 'r') as pdb_file:
        lines = pdb_file.readlines()

    molecules = []  # Lista per raccogliere le pose di ogni molecola
    current_molecule = []  # Buffer per le pose della molecola attuale
    pose_count = 0
    molecule_count = 0

    for line in lines:
        current_molecule.append(line)
        if line.strip() == "ENDMDL":
            pose_count += 1
            if pose_count % 9 == 0:  # Ogni 9 pose, salva la molecola
                molecules.append(current_molecule)
                current_molecule = []
                molecule_count += 1
                save_molecule(molecules, output_folder, molecule_count)
                molecules = []  # Reset della lista

    # Salva l'ultima molecola rimasta se non è stata salvata
    if current_molecule:
        molecules.append(current_molecule)
        molecule_count += 1
        save_molecule(molecules, output_folder, molecule_count)

def save_molecule(molecule_poses, output_folder, count):
    # Controlla se il file contiene solo "END" o "END\t"
    if len(molecule_poses) == 1 and (molecule_poses[0][0].strip() == "END" or molecule_poses[0][0].strip() == "END\t"):
        print(f"Il file molecule_{count}.pdb contiene solo 'END' o 'END\\t'. Non verrà salvato.")
        return

    filename = os.path.join(output_folder, f"molecule_{count}.pdb")
    with open(filename, 'w') as out_file:
        for pose in molecule_poses:
            out_file.writelines(pose)
    print(f"Salvato {filename}")

# Utilizzo
output_folder = "output_molecules"
split_pdb("docked.pdb", output_folder)