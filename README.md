# map_smallmol-FunFams
Piece of Python code for the code review meeting.

## What the code does
This code is adapted from the work by Felix Kruger at ChEMBL with mapping small molecules to Pfam domains. 
I've tweaked his code to make it work with FunFams instaed of Pfam domains.

## Requirements
* Access to CATH:FunFams Oracle database within Python (cx_Oracle). As far as I know cx_Oracle is available on bertha and bsmlx65

* Access to ChEMBL (MySQL): A local install of ChEMBL 19 is ready to use on bsmlx65

### Example of use
python masterFunFam.py

You should get the file 'map_funfam.txt' on the same directory you run the script.

## What I would like to get from this review
The code works fine but it does a few things I don't need, like checking for white-listed / black-listed Pfam domains.
It would be great to have some ideas on how to remove those. Otherwise, a double checking of errors (i.e. getting wrong FunFams...)
or any suggestion on code cleaning, improvement, etc are very welcome.

