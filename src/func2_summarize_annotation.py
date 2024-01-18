import json 
import os 

def func2_summarize_annotation(blastpOutput, outputPath):

    getHERVregion(blastpOutput, outputPath)

    HERVregion = f"{outputPath}/output_dict.json"

    createGTF(output_dict, outputPath)

def getHERVregion(blastpOutput, outputPath):

    output_dict = {}

    with open('../data/HERV_ORF_dict', 'r') as file:
        HERV_ORF_dict = json.load(file)

    with open(blastpOutput, 'r')as file:
        index = 0;
        for line in file:

            fields = line.strip().split('\t')
            key = index
            value = {
            "input_ID" : fields[0],
            "ORF_ID" : fields[1],
            "ORF" : fields[3],
            "len_alignment" : fields[4]
            }

            ## filter result 
            if float(fields[2]) != 100.000:
                continue
            elif fields[4] != fields[7]:
                continue

            else:
                index += 1
                value["HERV_ORF_dict"] = {
                "HERV region":HERV_ORF_dict[value["ORF_ID"]]["ID"],
                "ORF":HERV_ORF_dict[value["ORF_ID"]]["ORF"],
                "len_alignment":HERV_ORF_dict[value["ORF_ID"]]["len_alignment"],
                "strand":HERV_ORF_dict[value["ORF_ID"]]["strand"],
                "HERV_region_location":HERV_ORF_dict[value["ORF_ID"]]["region"]
                }

                region_start = int(HERV_ORF_dict[value["ORF_ID"]]["region"][0])
                region_end = int(HERV_ORF_dict[value["ORF_ID"]]["region"][1])
                ORF_start = int(fields[11])
                ORF_end = int(fields[12])
                aa_length = int(fields[4])

                start = region_start
                end = region_end

                ## if alignment length < ORF total length, need to recalculate the position of the region 
                if ( int(fields[4]) < int(fields[10]) ):
                    ## if start > end, it is reverse complement case
                    if region_start < region_end:
                        start = (ORF_start - 1) * 3 + region_start
                        end = start + aa_length * 3 - 1
                    else:
                        start = region_start - (ORF_start - 1) * 3
                        end = start - aa_length * 3 + 1

                value["location"] = [start, end]
                output_dict[index] = value

    with open(f"{outputPath}/output_dict.json", 'w') as file:
        json.dump(output_dict, file, indent=4)

def createGTF(output_dict, outputPath):
    
    with open(f"{outputPath}/quantification.gtf", 'w') as gtf_file:

        for i in range(1,len(output_dict) + 1):
            i = str(i)
            chromosome = output_dict[i]["ORF_ID"].split('_')[1]
            strand = '+'
            gene_id = output_dict[i]["ORF_ID"]
            start = output_dict[i]["location"][0]
            end = output_dict[i]["location"][1]

            ## determine start and end by reverse or forward strand
            if start > end:
                start = output_dict[i]["location"][1]
                end = output_dict[i]["location"][0]
                strand = '-'
            
            gtf_line = (
                f"{chromosome}\t"
                "blastp\t"
                "HERV\t"
                f"{start}\t"
                f"{end}\t"
                f".\t"
                f"{strand}\t"
                f".\t"
                f'gene_id "{i}_{output_dict[i]["input_ID"]}_{gene_id}"'
            )

            gtf_file.write(gtf_line + "\n")

    return None
