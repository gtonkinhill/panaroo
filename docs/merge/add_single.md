# Add a new genome to a Panaroo graph.

Having run Panaroo it is possible to add a single genome to the resulting pangenome graph. This can be useful when we are interested in the annotations of a new genome but do not wish to rerun the entire Panaroo algorithm. This added functionality is thanks to [Daniel Anderson](https://github.com/Danderson123).

As an example we initially run Panaroo on the species *Klebsiella quasipneumoniae* (KpII) dataset from [Holt et al., 2015](https://www.pnas.org/content/112/27/E3574) and query a single genome from the sister species *Klebsiella variicola* (KpIII).

```
panaroo -i kpII/*.gff -o kp_II_panaroo -t 24 --clean-mode sensitive
```

we can now add a single genome **5150_1#7.gff** to the graph by running

```
panaroo-integrate -d kp_II_panaroo/ -i kpIII/5150_1#7.gff -t 24 -o updated_output
```

The updated graph as well as updated gene presence absence files are available in the folder `updated_output`