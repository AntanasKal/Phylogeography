# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:46:53 2019

@author: Antanas
"""
#import dendropy



def write_phyrex_input(tree, i, input_string="output/phyrex_input/" , output_string="output/phyrex_output/"):
    write_xml(tree, i, input_string, output_string)
    write_phyrex_tree(tree, i, input_string, output_string)
    write_phyrex_nexus(tree, i, input_string, output_string)
    write_phyrex_coord(tree, i, input_string, output_string)
    
    
    
def write_xml(tree, i, input_string, output_string):
    file = open(input_string+"phyrex"+str(i)+".xml", "w")
    file.write("""<phyrex run.id=""" +str(i)+""" output.file="""+'"'+output_string+"out" +'"'+""" mcmc.chain.len="1E+6" mcmc.sample.every="1000"
        mcmc.print.every="1000" mcmc.burnin="10000" mutmap="no" ignore.sequences="yes">

  <!-- Tree topology -->
  <topology>
  """)
    file.write('\t <instance id="T1" init.tree="user" file.name="'+input_string+'phyrex_tree'+str(i)+'.txt" optimise.tree="no"/>\n')  
    file.write('</topology>\n')
    file.write("""
               <!-- Model of rate variation across lineages -->
  <lineagerates model="lognormal"/>

  <!-- Average (clock) rate of substitution -->
  <clockrate value="1E-4"/>
  
  
  <!-- Substitution model -->
  <ratematrices id="RM1">
    <instance id="M1" model="HKY85" optimise.tstv="no" tstv="4.0"/>
  </ratematrices>

  
    <!-- Freerate model of variation of rates across sites -->
  <siterates id="SR1">
    <instance id="R1" init.value="1.0"/>
    <weights  id="D1" family="gamma" optimise.freerates="no">
      <instance appliesto="R1" value="1.00"/>
    </weights>
  </siterates>

  <!-- Nucleotide frequencies -->
  <equfreqs id="EF1">
    <instance id="F1" optimise.freqs="no"/>
  </equfreqs>


  <!-- Vector of edge lengths -->
  <branchlengths id="BL1" >
    <instance id="L1" optimise.lens="no"/>
  </branchlengths>

  <!-- Model assembly -->
""")
    file.write('  <partitionelem id="partition1" file.name="'+input_string+'phyrex_nexus'+str(i)+'.nxs" data.type="nt" interleaved="no">\n')
    file.write("""    <mixtureelem list="T1"/>
    <mixtureelem list="M1"/>
    <mixtureelem list="F1"/>
    <mixtureelem list="R1"/>
    <mixtureelem list="L1"/>
  </partitionelem>
  """)
    
    file.write('  <coordinates id="coordinates" file.name="'+input_string+'phyrex_coord'+str(i)+'.txt"/>\n')
    index = 1
    for leaf in tree.leaf_node_iter():
        file.write('\t<clade id="clad'+str(index)+'">\n')
        file.write('\t\t<taxon value="'+leaf.taxon.label+'"/>\n')
        file.write('\t</clade> \n')
        
        file.write('\t<calibration id="cal'+str(index)+'">\n')
        file.write('\t\t<lower>'+str(leaf.time)+'</lower>\n')        
        file.write('\t\t<upper>'+str(leaf.time)+'</upper>\n')
        file.write('\t\t<appliesto clade.id="clad'+str(index)+'"/>\n')        
        file.write('\t</calibration> \n')
        index = index+1
    
    
    file.write('</phyrex>')
    file.close()
    
def write_phyrex_tree(tree, i, input_string, output_string):
    tree.write(path="output/phyrex_input/phyrex_tree"+str(i)+".txt", schema="newick", suppress_internal_taxon_labels=True)
    
def write_phyrex_coord(tree, i, input_string, output_string):
    file = open(input_string+"phyrex_coord"+str(i)+".txt", "w")
    file.write("# state.name lon lat\n")
    for leaf in tree.leaf_node_iter():
        file.write(leaf.taxon.label+' '+str(leaf.X)+' '+str(leaf.Y)+'\n')
    file.write("""|SouthWest| -50 -50
|NorthEast| 50 50 """)  
    file.close()
    
def write_phyrex_nexus(tree, i, input_string, output_string):
    num_leaves = 0
    for leaf in tree.leaf_node_iter():
        num_leaves=num_leaves+1
    
    file = open(input_string+"phyrex_nexus"+str(i)+".nxs", "w")
    
    file.write('#NEXUS\n')
    file.write('BEGIN DATA:\n')
    file.write('\tDIMENSIONS NTAX='+str(num_leaves)+' NCHAR='+str(12)+';\n')
    file.write(""" FORMAT DATATYPE=DNA INTERLEAVE;
MATRIX
""")
    for leaf in tree.leaf_node_iter():
        file.write(leaf.taxon.label+" CCAAAAGATAAT\n")
        
    file.write("""
;
END;""")
    file.close()