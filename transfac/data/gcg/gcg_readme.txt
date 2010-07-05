GCG-Files of TRANSFAC Professional 9.2 - (Releasedate 2005-06-30)

File Name       Species/Group
-------------------------------------------------------------------------
tfsites.arts	Artificial sequences
tfsites.cons	Consensus sequences
tfsites.fungi	Fungi
tfsites.insect  Arthropoda (Insects)
tfsites.other	Other (NOT arts, cons, fungi, insect, plant, vert, viral)
tfsites.plant	Plantae
tfsites.vert	Vertebrata
tfsites.viral	Viridae


How to use the TRANSFAC GCG formated files with e.g. program 
"findpatterns".

Unpack the files in the GenRunData directory in GCG. GCG has an older 
version available as "tfsites.dat" but in a different directory so it is 
harder to use.

The following text is put at the top of each file. The only essential
part is the line with ".." anything above that is free text, except
that it must never contain the string ".."

Example:

To use the vertebrate-file in findpatterns as the pattern file:

      % findpatterns -data=tfsites.vert


==========================start of text ======================================
TRANSFAC Release 9.2
GCG format sites files for artificial sequences, consensus sequences, fungi,
insects, other, plants, vertebrates, viruses.

TRANSFAC GCG data files are:

tfsites.arts
tfsites.cons
tfsites.fungi
tfsites.insect
tfsites.other
tfsites.plant
tfsites.vert
tfsites.viral

To use this file in findpatterns as the pattern file:

      % findpatterns -data=tfsites.vert

..
==========================  end of text ======================================

If you have problems using the files, please E-Mail to:

E-Mail:         support@biobase.de
URL:            http://www.biobase.de/

Mail:           BIOBASE GmbH
                Halchtersche Strasse 33
                38304 Wolfenbuettel
                Germany

Tel.:           +49 (0)5331 8584-0
Fax:            +49 (0)5331 8585-70

IRE/2005-06-24
