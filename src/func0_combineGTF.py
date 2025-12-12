import pandas as pd
import sys

def func0_combineGTF(gtfPath, outputPath):
    """
	Combine gtf files
	Input : 
		gtfPath : path to the csv file with paths of the gtf files
	Output : a combined gtf file 
	"""
    try:
        output_gtf = f"{outputPath}/combined_quantification.gtf"
        gtfPath_df = pd.read_csv(gtfPath, header=None)
        if gtfPath_df.shape[1] != 1:
            raise ValueError("Check the format of the input csv file")
        column_name = gtfPath_df.columns[0]
        column_data = gtfPath_df[column_name]
        gft_list = column_data.tolist()

        region_list = []

        with open(output_gtf, 'w') as f:

            for gtf_file in gft_list:
                with open(gtf_file, 'r') as gtf_f:
                    full_line = [line.strip() for line in gtf_f]
                    
                annot_line = [line.split('\t') for line in full_line]

                for i in range(len(full_line)):
                    chr = annot_line[i][0]
                    start = annot_line[i][3]
                    end = annot_line[i][4]
                    strand = annot_line[i][6]
                    region = f"{chr}_{start}_{end}_{strand}"
                    if region not in region_list:

                        region_list.append(region)
                        f.write(f"{full_line[i]}\n")

        return None
    
    except FileNotFoundError as err:
        print("Error: The file was not found. Please check the file path.")
        return err
    except ValueError as err:
        print(f"Error: please check your csv input: {err}", file=sys.stderr)
        return err
    except Exception as err:
        print(f"Error: {err}", file=sys.stderr)
        return err

