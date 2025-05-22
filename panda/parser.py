def parse_C6_C12(text: str, atom_types):
    lines = text.split("\n")
    flag = False

    params = dict(zip(atom_types, [[-1, -1] for _ in atom_types]))
    for line in lines:
        if line == "[ atomtypes ]":
            flag = True
            continue
        elif line != "":
            if line[0] == ";":
                continue
        elif line == "" and flag:
            break

        if flag:
            context = line.split()

            assert len(context) >= 5, f"Wrong line format: {line}"

            if context[0] in atom_types:
                params[context[0]][0] = float(context[4].strip(";"))
                params[context[0]][1] = float(context[5].strip(";"))

    # Checking case when atom_type is not founded in topol
    for key, val in params.items():
        assert not (val[0] == -1 or val[1] == -1), f"{key} is not founded in the text"

    return params
