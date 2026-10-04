"""Generate the standalone PyMOL runtime script from a prepared scene."""
import textwrap

_PYMOL_SCENE = r'''
from pathlib import Path
import json
import textwrap
import traceback
import numpy as np
from pymol import cmd

def chain_selection(name, author, label):
    selection = f"(model {name} and chain {json.dumps(author)} and polymer)"
    labelled = selection + f" and segi {json.dumps(label)}"
    if cmd.count_atoms(labelled):
        selection = labelled
    if not cmd.count_atoms(selection):
        raise ValueError(f"PyMOL author chain selection has no polymer atoms: {name}, {author}, label {label}")
    return selection

try:
    cmd.reinitialize()
    cmd.set("cif_use_auth", 1)
    cmd.set("matrix_mode", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("opaque_background", 1)
    cmd.bg_color("white")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", .16)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_size", 24)
    cmd.set("label_font_id", 7)
    cmd.set("float_labels", 1)
    cmd.set("label_distance_digits", 3)
    cmd.viewport(scene["width"], scene["height"])
    focus = []
    highlights = []
    atom_counts = {}
    highlighted_coords = []
    transformed_objects = []
    for molecule in scene["molecules"]:
        name = molecule["name"]
        cmd.load(molecule["path"], "raw_input", multiplex=0)
        cmd.create(name, "raw_input", source_state=molecule["state"], target_state=1)
        cmd.delete("raw_input")
        atom_counts[name] = cmd.count_atoms(name)
        if not atom_counts[name]:
            raise ValueError(f"No atoms loaded from {molecule['path']} state {molecule['state']}")
        if molecule["transform"] is not None:
            before = cmd.get_coords(name, state=1)
            cmd.transform_object(name, molecule["transform"], state=1, homogenous=1, transpose=0)
            after = cmd.get_coords(name, state=1)
            matrix = np.array(molecule["transform"]).reshape(4, 4)
            expected = before @ matrix[:3, :3].T + matrix[:3, 3]
            if not np.allclose(after, expected, atol=2e-4, rtol=1e-6):
                raise ValueError(f"PyMOL whole-molecule transform disagrees with saved rigid transform: {molecule['chain_key']}")
            transformed_objects.append(name)
        cmd.hide("everything", name)
        cmd.color("grey70", name)
        cmd.show("cartoon", f"model {name} and polymer")
        if scene["overlay"]:
            cmd.set("cartoon_transparency", .8, name)
        elif scene["name"] == "7baa_context":
            cmd.set("cartoon_transparency", .35, name)
        selected = chain_selection(name, molecule["auth_asym_id"], molecule["label_asym_id"])
        if molecule.get("domain_resi"):
            cmd.hide("everything", name)
            selected += " and resi " + json.dumps(molecule["domain_resi"])
            cmd.show("cartoon", selected)
        focus.append(selected)
        cmd.color(molecule["colour"], selected)
        cmd.set("cartoon_transparency", 0, selected)
        for author in molecule.get("partner_chains", []):
            partner = f"model {name} and chain {json.dumps(author)} and polymer"
            if not cmd.count_atoms(partner):
                raise ValueError(f"Requested partner chain absent: {name}, {author}")
            cmd.show("cartoon", partner)
            cmd.color("orange", partner)
        for peptide in ([] if molecule.get("domain_resi") else molecule["peptides"]):
            peptide_selection = chain_selection(name, peptide["auth_asym_id"], peptide["label_asym_id"])
            cmd.color("orange", peptide_selection)
            cmd.show("sticks", peptide_selection)
        if not molecule.get("domain_resi"):
            cmd.show("sticks", f"model {name} and organic")
        cmd.color("grey50", f"model {name} and organic and elem C")
        cmd.color("red", f"model {name} and organic and elem O")
        cmd.color("blue", f"model {name} and organic and elem N")
        for code in molecule["ligand_codes"]:
            ligand_selection = f"model {name} and organic and resn {json.dumps(code)}"
            cmd.show("sticks", ligand_selection)
            if scene["name"] == "reference_context":
                cmd.color("forest", ligand_selection + " and elem C")
                cmd.set("stick_radius", .22, ligand_selection)
            ligand_atoms = cmd.index(ligand_selection)
            if ligand_atoms:
                object_name, atom_index = ligand_atoms[0]
                cmd.label(f"model {object_name} and index {atom_index}", repr(code))
        if not molecule.get("domain_resi"):
            cmd.show("sticks", f"model {name} and polymer.nucleic")
        if molecule["highlight_resi"] is not None:
            highlight = selected + f" and resi {json.dumps(molecule['highlight_resi'])} and name CA and alt {json.dumps(molecule['highlight_altloc'])}"
            # The saved alternate choice is the one represented in the mapping.
            # More than one CA would make a distance label ambiguous.
            if cmd.count_atoms(highlight) != 1:
                raise ValueError(f"Outlier CA selection must contain exactly one atom: {highlight}")
            cmd.show("spheres", highlight)
            cmd.set("sphere_scale", 1, highlight)
            cmd.label(highlight, repr(molecule["highlight_label"]))
            cmd.set("label_position", [0, 3, 0], highlight)
            highlighted_coords.append(cmd.get_coords(highlight, state=1)[0].tolist())
            highlights.append(highlight)
    measured_displacement = None
    if len(highlights) == 2:
        measured_displacement = float(np.linalg.norm(np.array(highlighted_coords[1]) - highlighted_coords[0]))
        if not np.isclose(measured_displacement, scene["displacement_A"], atol=2e-4, rtol=1e-6):
            raise ValueError(f"Rendered CA distance {measured_displacement} differs from saved displacement {scene['displacement_A']}")
        cmd.distance("saved_local_displacement", highlights[0], highlights[1], cutoff=measured_displacement + 1)
        cmd.set("dash_color", "black", "saved_local_displacement")
        cmd.set("label_position", [0, -3, 0], "saved_local_displacement")
    orientation = " or ".join(focus) if scene["overlay"] else "all"
    if highlights:
        # Orient by the measured pair so its separation lies in the image
        # plane rather than disappearing behind another helix along view depth.
        orientation = " or ".join(highlights)
    cmd.orient(orientation)
    cmd.zoom("all", buffer=5, complete=1)
    if scene.get("view") is not None:
        cmd.set_view(scene["view"])
    cmd.deselect()
    cmd.png(scene["png"], width=scene["width"], height=scene["height"], dpi=180, ray=1, quiet=1)
    if not Path(scene["png"]).exists():
        raise RuntimeError("PyMOL returned without writing the requested PNG")
    # Keep the actual molecular render intact and append an explanatory banner.
    # This annotation code lives in the PML as well as in the generated PNG.
    from PIL import Image, ImageDraw, ImageFont
    image = Image.open(scene["png"]).convert("RGB")
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 19)
    except OSError:
        font = ImageFont.load_default(size=19)
    lines = textwrap.wrap(scene["caption"], width=max(45, image.width // 11))
    banner_height = 35 + 27 * len(lines)
    annotated = Image.new("RGB", (image.width, image.height + banner_height), "white")
    annotated.paste(image, (0, 0))
    draw = ImageDraw.Draw(annotated)
    if scene.get("legend"):
        legend = scene["legend"]
        legend_width = int(max(draw.textlength(line, font=font) for line in legend)) + 30
        draw.rectangle((10, 10, 10 + legend_width, 25 + 27 * len(legend)), fill="white")
        for index, line in enumerate(legend):
            colour = "#228833" if line.startswith("Green") else "black"
            draw.text((20, 20 + 27 * index), line, fill=colour, font=font)
    if measured_displacement is not None:
        # Ray-traced helices can occlude native labels. Keep the numerical
        # result and both highlighted identities readable in the empty corner.
        callout = [molecule["highlight_label"] for molecule in scene["molecules"]]
        callout.append(f"CA displacement = {measured_displacement:.3f} A")
        callout_width = int(max(draw.textlength(line, font=font) for line in callout)) + 30
        draw.rectangle((10, 10, 10 + callout_width, 25 + 27 * len(callout)), fill="white")
        for index, line in enumerate(callout):
            colour = "#0077bb" if index == 0 else "#aa0099" if index == 1 else "black"
            draw.text((20, 20 + 27 * index), line, fill=colour, font=font)
    for index, line in enumerate(lines):
        draw.text((20, image.height + 15 + 27 * index), line, fill="black", font=font)
    annotated.save(scene["png"])
    Path(scene["details"]).write_text(json.dumps({"atom_counts": atom_counts, "pymol_version": cmd.get_version()[0], "whole_molecule_transform_checked": bool(transformed_objects), "transformed_objects": transformed_objects, "highlighted_CA_coords_A": highlighted_coords, "measured_CA_displacement_A": measured_displacement}, indent=2) + "\n")
except Exception:
    Path(scene["error"]).write_text(traceback.format_exc())
    cmd.quit(1)
else:
    cmd.quit(0)
'''


def write_scene_script(scene, path):
    script = "# Reproduce with: pymol -cq " + str(path) + "\npython\n"
    script += "scene = " + repr(scene) + "\n" + textwrap.dedent(_PYMOL_SCENE) + "\npython end\n"
    path.write_text(script, encoding="utf-8")
