import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("--WIDTH_X", type=float, required=True)
parser.add_argument("--WIDTH_Y", type=float, required=True)
parser.add_argument("--H", type=float, required=True)
parser.add_argument("--phi", type=float, required=True)
parser.add_argument("--delta_list", type=str, required=True)
parser.add_argument("--interface_type", type=str, required=True)
parser.add_argument("--extention", type=str, required=True)
parser.add_argument("--folder", type=str, required=True)
parser.add_argument("--iterations", type=int, required=True)

args = parser.parse_args()
args_vars = vars(args)
args_vars["l"] = args.WIDTH_X / args.H

with open("scripts/run_delta_template.sh", "r") as template_file:
    template_content = template_file.read()

script_content = template_content.format(**args_vars)

# Define a unique filename for each script
script_filename = os.path.join("scripts", f"run_delta_{args.interface_type}.sh")

# Write the generated content to the new script file
with open(script_filename, "w") as script_file:
    script_file.write(script_content)
