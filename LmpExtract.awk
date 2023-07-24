#!/bin/awk -f
# USAGE: Execute this script with the 'output' file as an argument
# EXAMPLE: ./LmpExtract output


# Search every ocurrence and count it
/^Step/ {
  count++
}

# If is the fifth ocurrence and the line begins with a space save it
count==5 && /^\s+/ {

  # Remove the first spaces (this is for better formating)
  gsub (/^\s+/, "")

  # Save the data in this file
  print $0 > "output.dat"
}

