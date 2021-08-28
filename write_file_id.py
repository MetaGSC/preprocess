import csv

from helpers import timestamp

def write_file_id(batch_i, input_arr, target_path, err_file):
    try:
        
        with open(target_path, "a", encoding='UTF8') as f:
            writer = csv.writer(f)
            for n, frag in input_arr:
                data = [batch_i, n, frag["id"].strip().split()[0]]
                if("name" in frag.keys()):
                    data.append(frag["name"])
                writer.writerow(data)
    except Exception as err:
        with open(err_file, 'a') as fout:
            fout.write(f"{timestamp()} Error writing file id into csv for the batch:{batch_i} {err}\n")
