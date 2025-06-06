import os
import sys
import assists

def gen_table_summary(full_info):
    # Encode all the set images
    dirpath = os.path.dirname(os.path.abspath(full_info))
    workdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    template_html = workdir + "/assets/table_template.html"
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

    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        print(f"Replacing {key} to {value} in the HTML")
        template = template.replace(key, str(value))

    output_file = dirpath + "/table_summary.html"
    with open(output_file, "w") as outfile:
        outfile.write(template)

"""
<button class="accordion">py_species1_ph</button>
    <div class="panel">
        <p>
            py_species1status_ph
        </p>
        <p style="padding-bottom: 10px;">
            py_species1zscore_ph
        </p>
    </div>
"""

if __name__ == "__main__":
    gen_table_summary(sys.argv[1])
