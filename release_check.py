# check that release version across CHANGELOG.md, Project.toml, and git tag all agree
# if executed as `release_check.py publish [remote name]` will push the latest git tag
# to the specified remote

import os
import subprocess
import sys

def get_release_projecttoml():
    with open("Project.toml", "r") as f:
        for line in f:
            if line.strip().lower().startswith("version"):
                lsp = line.split("=")
                v = lsp[1].replace('"', '')
                return [int(x) for x in v.split('.')]

def get_release_changelog():
    with open("CHANGELOG.md", "r") as f:
        for line in f:
            if line.strip().lower().startswith("## v"):
                lsp = line.split()
                v = lsp[1].replace('v', '')
                return [int(x) for x in v.split('.')]
    # Changelog

    ## v0.6.2
def get_release_gittag():
    output = subprocess.run(["git", "tag"], stdout=subprocess.PIPE)
    decoded = output.stdout.decode("utf-8")
    tags = [t.strip().lstrip('v') for t in decoded.strip().split('\n') if t.strip()]
    
    # Filter out tags that can't be split into version numbers
    valid_tags = []
    for t in tags:
        try:
            parts = [int(x) for x in t.split('.')]
            valid_tags.append((t, parts))
        except ValueError:
            continue  # Skip malformed tags

    if not valid_tags:
        raise ValueError("No valid version tags found")

    # Sort by version number
    valid_tags.sort(key=lambda x: x[1])
    return [int(x) for x in valid_tags[-1][0].split('.')]

def find_docs_todos():
    todos = []
    for root, _, files in os.walk(os.path.join("docs", "src")):
        for fname in files:
            if not fname.endswith(".md"):
                continue
            fpath = os.path.join(root, fname)
            with open(fpath, "r") as f:
                for lineno, line in enumerate(f, start=1):
                    if "TODO" in line:
                        todos.append((fpath, lineno, line.strip()))
    return todos

v_projecttoml = get_release_projecttoml()
v_changelog = get_release_changelog()
v_gittag = get_release_gittag()


str_v_pt = '.'.join([str(x) for x in v_projecttoml])
str_v_cl = '.'.join([str(x) for x in v_changelog])
str_v_gt = '.'.join([str(x) for x in v_gittag])

agree = True
if v_projecttoml != v_changelog:
    agree = False
    print(f"Version numbers in Project.toml ({str_v_pt}) and CHANGELOG.md ({str_v_cl}) are different!!!")

if v_gittag != v_projecttoml:
    agree = False
    print(f"Version numbers in last git tag ({str_v_gt}) and Project.toml ({str_v_pt}) are different!!!")

if v_gittag != v_changelog:
    agree = False
    print(f"Version numbers in last git tag ({str_v_gt}) and CHANGELOG.md ({str_v_cl}) are different!!!")

docs_todos = find_docs_todos()
if docs_todos:
    agree = False
    print("Found hanging TODOs in docs/src:")
    for fpath, lineno, line in docs_todos:
        print(f"  {fpath}:{lineno}: {line}")

publish = False
remote_name = "unknown"
if len(sys.argv) > 2:
    if sys.argv[1].lower() == "publish":
        publish = True
        remote_name = sys.argv[2]

if not agree:
    print("Versions disagree!!!")
else:
    print("OK")
    if publish:
        str_v_gt = "v" + str_v_gt
        y_n = input(f"Pushing {str_v_gt} to {remote_name}, Y to continue")
        if y_n.strip().lower() == "y":
            output = subprocess.run(["git", "push", remote_name, str_v_gt])
            output.check_returncode()