import csv

from helpers import timestamp

def write_file_id(input, target_path, err_file):
    try:
        n, frag = input
        data = [n, frag["id"].strip().split()[0]]
        if("name" in frag.keys()):
            data.append(frag["name"])

        with open(target_path, "a", encoding='UTF8') as f:
            writer = csv.writer(f)
            writer.writerow(data)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing file id into csv for the file:{n} {err}\n")