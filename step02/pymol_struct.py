import os
import pymol

# launching Pymol
pymol.finish_launching()

# set path
parent_folder = "~/output/"
pdb_files = [f for f in os.listdir(parent_folder) if f.endswith(".pdb")]


for pdb_file in pdb_files:

    # check rank1 result
    if "rank_001" in pdb_file:

        # load pdb
        pdb_path = os.path.join(parent_folder, pdb_file)
        pymol.cmd.load(pdb_path)
        
        # style for molecules visualization
        pymol.cmd.show("cartoon")
        pymol.cmd.color("palecyan")
        pymol.cmd.color("red", "ss h")
        pymol.cmd.color("gold", "ss s")
        
        # Render image
        pymol.cmd.ray(800, 600)
        
        # filename
        file_prefix = '_'.join(pdb_file.split('_')[:2])
        
        # Build a path to save the image
        image_name = f"{file_prefix}_st.png"
        image_path = os.path.join("~/HMPA_resource/images/pdb", image_name)
        os.makedirs(os.path.dirname(image_path), exist_ok=True)
        
        # save
        pymol.cmd.png(image_path)
        pymol.cmd.delete("all")

# close pymol
pymol.cmd.quit()