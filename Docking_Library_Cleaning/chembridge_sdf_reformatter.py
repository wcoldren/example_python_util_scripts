#!/usr/bin/env python

def reformat_chembridge_sdf(chembridge_sdf: str, output_sdf: str):
    """Take a ChemBridge SDF from their website and reformat it for use in
    both RDKit and Schrodinger


    Keyword arguments:
    chembridge_sdf -- type: str -- filename of the sdf input
    output_sdf -- type: str -- output sdf filename
    """
    raw_sdf = open(chembridge_sdf).read().splitlines()

    for idx, line in enumerate(raw_sdf):
        # Remove ISIS lines, they serve no purpose, also use them to anchor
        # The top line of the SDF description which will be the compound ID.
        if "ISIS" in line:
            header_line = idx - 1
            raw_sdf[idx] = ''
        # Grab the ChemBridge ID number and substitute it into the appropriate
        # spot in the list.
        elif "<ID>" in line:
            chembridge_id = raw_sdf[idx + 1]
            raw_sdf[header_line] = chembridge_id

    with open(output_sdf, 'w') as f:
        for line in raw_sdf:
            # If you are confused look up f-strings, new since python 3.6.
            f.write(f"{line}\n")


def main():
    reformat_chembridge_sdf("./Input/exp_sdf.sdf", "./Sample_Out/exp_sdf_reformatted.sdf")


if __name__ == "__main__":
    main()
