import sys

start_key = "$SERIALBOX START FV NML"
stop_key = "$SERIALBOX STOP FV NML"


def main(log_file: str):
    with open(log_file, "r") as f:
        lines = f.readlines()

    line_start = -1
    for i, line in enumerate(lines):
        if start_key in line:
            line_start = i

    with open("input.nml", "w") as nml_f:
        nml_f.write("&fv_core_nml\n")

        for line_index in range(line_start, len(lines)):
            line = lines[line_index]
            if stop_key in line:
                break
            line = line.replace("FV3", "").replace("\n", "")
            elts = line.split(":")
            if len(elts) == 2:
                if elts[1] != "":
                    value = elts[1].strip()
                    if value == "F":
                        value = ".false."
                    elif value == "T":
                        value = ".true."
                    key = elts[0].strip()
                    nml_f.write(f" {key} = {value}\n")
        nml_f.write(" layout = 1,1\n")
        nml_f.write("/\n")


if __name__ == "__main__":
    main(sys.argv[1])
