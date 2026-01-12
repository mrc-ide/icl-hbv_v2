## Script to be run in src/raw_params. Takes the "params.txt" file and pulls out the country name and ISO

params_file = open("params.txt","r")
##params_file = open("temp1.txt","r")
## Creates an array where each entry is a country:
params_raw = params_file.read().rstrip().split(" COUNTRY: ")[1:]

outstring = ""
for country_data in params_raw:
    split_data = country_data.splitlines()
    [ISO_data, country_data] = split_data[:2]
    country_data = country_data.rstrip().replace('"country_name": "','').replace('",','')
    outstring += ISO_data[0:3] + ","+ country_data +"\n"
outstring = outstring.rstrip()

outfile = open("processed_params/ListOfCountryNamesAndISOS.csv","w")
outfile.write(outstring)
outfile.close()
