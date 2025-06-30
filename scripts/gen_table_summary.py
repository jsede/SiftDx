import os
import sys
import assists
import pandas as pd

pathogen_db = pd.read_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "databases/pathogen_list.csv"), header = 0)
contam_db = pd.read_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "databases/known_contaminants.csv"), header = 0)
parasite_list = pd.read_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "databases/parasites.txt"), sep="\t",header = 0)

def gen_table_summary(full_info):
    # file set-up    
    dirpath = os.path.dirname(os.path.abspath(full_info))
    workdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    template_html = workdir + "/assets/table_template.html"
    
    # Encode all the set images
    img_names = [
        "bacteria.png", 
        "virus.png", 
        "fungi.png", 
        "parasite.png"
        ]
    
    base64_images = []

    for img in img_names:
        img_path = os.path.join(workdir, "assets", img)
        base64_image = assists.base64encode(img_path)
        variable_name = "base64_" + os.path.splitext(img)[0]
        locals()[variable_name] = base64_image
        base64_images.append(variable_name)

    replace_dict = {
        "py_bacteria_icon_ph": locals()[base64_images[0]],
        "py_virus_icon_ph": locals()[base64_images[1]],
        "py_fungi_icon_ph": locals()[base64_images[2]],
        "py_parasite_icon_ph": locals()[base64_images[3]],
    }
    # load the sample data
    zscore_df = pd.read_csv(full_info, sep="\t", header=0)

    # apply all the labels.
    zscore_df[['taxon_label', "taxon_status"]] = zscore_df.apply(
        lambda row: assists.flag_pathogens(row, pathogen_db, contam_db),
        axis=1,
        result_type='expand'
    )
    # grab bacterial info
    bacteria_data = assists.taxon_split(zscore_df, "Bacteria", "superkingdom")
    bacteria_html = assists.html_loop(bacteria_data)
    replace_dict["py_bacteria_ph"] = bacteria_html

    # grab viral info
    viral_data = assists.taxon_split(zscore_df, "Viruses", "superkingdom")
    viral_html = assists.html_loop(viral_data)
    replace_dict["py_virus_ph"] = viral_html
    
    # grab fungal info
    fungal_data = assists.taxon_split(zscore_df, "Fungi", "kingdom")
    fungal_html = assists.html_loop(fungal_data)
    replace_dict["py_fungi_ph"] = fungal_html

    # grab parasite info
    parasite_data = assists.get_parasites(zscore_df, parasite_list)
    parasite_html = assists.html_loop(parasite_data)
    replace_dict["py_parasite_ph"] = parasite_html

    # update the html with new rows for each pathogen etc.
    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        print(f"Replacing {key} to {value} in the HTML")
        if isinstance(value, list):  # for HTML fragments like pathogen rows
            template = template.replace(key, '\n'.join(value))
        else:  # for base64 strings
            template = template.replace(key, value)

    output_file = dirpath + "/table_summary.html"
    with open(output_file, "w") as outfile:
        outfile.write(template)



if __name__ == "__main__":
    gen_table_summary(sys.argv[1])
