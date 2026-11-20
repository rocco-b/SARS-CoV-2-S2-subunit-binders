input_file = "log.txt"
output_file = "out.txt"
delimiter = "-----+------------+------------+------------+----------"

count = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    lines = infile.readlines()
    
    for i in range(len(lines)):
        if lines[i].strip() == delimiter and i + 1 < len(lines):
            outfile.write(lines[i + 1])
            count += 1

print(f"Scritte {count} righe nel file di output")
