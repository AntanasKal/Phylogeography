# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:37:34 2019

@author: Antanas
"""

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max, help='sum the integers (default: find the max)')

args = parser.parse_args()
print(args.accumulate(args.integers))


import dendropy
import random
from dendropy.simulate import treesim

def simulate_brownian(t, var, dimension):
    for node in t.preorder_node_iter():
        if node.parent_node is None:
            node.X = float(0)
            node.displacementx = float(0)
            if dimension==2:
                node.Y = float (0)
                node.displacementy = float(0)
        else:
            node.displacementx = random.gauss(0, var*node.edge.length)
            node.X = node.parent_node.X+node.displacementx            
            if dimension==2:
                node.displacementy = random.gauss(0, var*node.edge.length)
                node.Y = node.parent_node.Y+node.displacementy            
    return t

def cal_times(t):
    for node in t.preorder_node_iter():
        if node.parent_node is None:
            node.time = 0
        else:
            node.time = node.parent_node.time+node.edge.length            
    return t    


def generate_tree(br, dr, num_extinct):
    t = treesim.birth_death_tree(birth_rate=br, death_rate=dr, num_extinct_tips=num_extinct, is_retain_extinct_tips=True, is_add_extinct_attr=True)
    #t.print_plot()    
    
    index = 0
    namespace = [];
    
    for node in t.preorder_node_iter():
        index=index+1
        namespace.append("T"+str(index))
    
    #name all nodes instead of just leaves
    taxon_namespace = dendropy.TaxonNamespace(namespace)
    t.taxon_namespace=taxon_namespace
    index=0
    for node in t.preorder_node_iter():
        index=index+1
        node.taxon=t.taxon_namespace.get_taxon("T"+str(index))
    
    t =prune_nodes(t)
    
    #distance to root
    t=cal_times(t)
        
    return t

def prune_nodes(t):
    for leaf in t.leaf_node_iter():    
        if hasattr(leaf, 'is_extinct'):
            leaf.extinct_ancestor = True
        else:
            leaf.extinct_ancestor = False
        
    for node in t.postorder_node_iter(): 
        if not hasattr(node, 'extinct_ancestor'):
            child_extinct = False
            for child in node.child_node_iter():
                if child.extinct_ancestor:
                    child_extinct =True
            node.extinct_ancestor = child_extinct
    labels = set([taxon.label for taxon in t.taxon_namespace
        if not t.find_node_for_taxon(taxon).extinct_ancestor])
    t1 = t.extract_tree_without_taxa_labels(labels=labels)

    return t1




def generate_coalescent_tree():
    num_tips = 20
    names = []
    for i in range(2*num_tips-1):
        names.append("T"+str(i))
    
#     print(names)
    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    current_nodes = []
    for i in range(num_tips):
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(i)))
        current_nodes.append(node)
        node.age = 0
        
    
    
    for merges in range(num_tips-1):
        time_to_coalescent=random.expovariate(len(current_nodes)*(len(current_nodes)-1)/2)
        time_from_present=time_from_present+time_to_coalescent
        merging_branches = random.sample(range(len(current_nodes)),2)
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(merges+num_tips)))
        if merges == num_tips-2:
            node=tree.seed_node
            node.taxon=taxon_namespace.get_taxon("T"+str(merges+num_tips))
        node.age = time_from_present
        current_nodes[merging_branches[0]].edge.length=time_from_present-current_nodes[merging_branches[0]].age
        current_nodes[merging_branches[1]].edge.length=time_from_present-current_nodes[merging_branches[1]].age
        node.set_child_nodes([current_nodes[merging_branches[0]], current_nodes[merging_branches[1]]])
        
        current_nodes.pop(max(merging_branches))
        current_nodes.pop(min(merging_branches))
        current_nodes.append(node)
#     print(tree.as_string("newick"))
#     print(tree.as_ascii_plot(show_internal_node_labels=True, plot_metric='length'))
    tree=cal_times(tree)
    return tree



def generate_coalescent_nonultrametric_tree():
    lamb=1
    period_length=1.3
    num_tips_per_period = 5
    num_periods = 4
    num_tips = num_tips_per_period*num_periods
    names = []
    for i in range(2*num_tips-1):
        names.append("T"+str(i))
#     print(names)

    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    current_nodes = []
    index = 0
    
    for current_period in range(num_periods):
        time_from_present=current_period*period_length
        for i in range(num_tips_per_period):
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(index)))
            current_nodes.append(node)
            index= index+1
            node.age = time_from_present
        
        
        current_num_tips = len(current_nodes)
        
        for merges in range(current_num_tips-1):
            time_to_coalescent=random.expovariate(lamb*len(current_nodes)*(len(current_nodes)-1)/2)
            print(len(current_nodes))
            time_from_present=time_from_present+time_to_coalescent
            if current_period < num_periods-1 and time_from_present > (current_period+1)*period_length:
                break
            else:
                merging_branches = random.sample(range(len(current_nodes)),2)
                
                if merges == current_num_tips-2 and current_period==num_periods-1:
                    node=tree.seed_node
                    node.taxon=taxon_namespace.get_taxon("T"+str(index))
                else:
                    node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(index)))
                index=index+1
                    
                node.age = time_from_present
                current_nodes[merging_branches[0]].edge.length=time_from_present-current_nodes[merging_branches[0]].age
                current_nodes[merging_branches[1]].edge.length=time_from_present-current_nodes[merging_branches[1]].age
                node.set_child_nodes([current_nodes[merging_branches[0]], current_nodes[merging_branches[1]]])
        
                current_nodes.pop(max(merging_branches))
                current_nodes.pop(min(merging_branches))
                current_nodes.append(node)
    print(tree.as_string("newick"))
    print(tree.as_ascii_plot(show_internal_node_labels=True, plot_metric='length'))
    tree=cal_times(tree)
    return tree

tree=generate_coalescent_nonultrametric_tree()

for node in tree.preorder_node_iter():
    print("%s : %s : %s" % (node.taxon.label, node.time, node.age))


















def write_BEAST_xml(t, d, i, dimension):
    if dimension==2:
        write_BEAST_xml_dim_2(t, d, i)
    else:
        write_BEAST_xml_dim_1(t, d, i)
        
def write_BEAST_xml_dim_2(t, d, i):
    file = open("output8/beast"+str(i)+".xml","w")
    file.write('<?xml version="1.0" standalone="yes"?>\n')
    file.write('<beast version="1.10.4">\n')
    file.write('\t<taxa id="taxa">\n')
    for tax in d:
        file.write('\t\t<taxon id="'+tax.label+'">\n')
        file.write('\t\t\t<date value="'+str(t.find_node_for_taxon(tax).time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).X)+'\n')
        file.write('\t\t\t</attr>\n')
        file.write('\t\t\t<attr name="Y">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).Y)+'\n')
        file.write('\t\t\t</attr>\n')
        
        ##perhaps not needed?
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).X)+'\n')
        file.write('\t\t\t</attr>\n')
        file.write('\t\t\t<attr name="Y">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).Y)+'\n')
        file.write('\t\t\t</attr>\n')

        file.write('\t\t</taxon>\n')   
    file.write('\t</taxa>\n')  
    
    file.write('\t<newick id="startingTree">\n')
    file.write('\t\t'+t.as_string(schema="newick",suppress_rooting=True)+'\n')
    
    file.write('\t</newick>\n')
    

    
    file.write("""	<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>\n""")
    
    file.write("""	<!-- Statistic for sum of the branch lengths of the tree (tree length)       -->
	<treeLengthStatistic id="treeLength">
		<treeModel idref="treeModel"/>
	</treeLengthStatistic>

	<!-- Statistic for time of most recent common ancestor of tree               -->
	<tmrcaStatistic id="age(root)" absolute="true">
		<treeModel idref="treeModel"/>
	</tmrcaStatistic>

<!-- START Multivariate diffusion model                                      -->

	<multivariateDiffusionModel id="X.diffusionModel">
		<precisionMatrix>
			<matrixParameter id="X.precision">
				<parameter id="X.precision.col1" value="0.05"/>
			</matrixParameter>
		</precisionMatrix>
	</multivariateDiffusionModel>

	<multivariateWishartPrior id="X.precisionPrior" df="1">
		<scaleMatrix>
			<matrixParameter>
				<parameter value="1.0"/>
			</matrixParameter>
		</scaleMatrix>
		<data>
			<parameter idref="X.precision"/>
		</data>
	</multivariateWishartPrior>

	<multivariateDiffusionModel id="Y.diffusionModel">
		<precisionMatrix>
			<matrixParameter id="Y.precision">
				<parameter id="Y.precision.col1" value="0.05"/>
			</matrixParameter>
		</precisionMatrix>
	</multivariateDiffusionModel>

	<multivariateWishartPrior id="Y.precisionPrior" df="1">
		<scaleMatrix>
			<matrixParameter>
				<parameter value="1.0"/>
			</matrixParameter>
		</scaleMatrix>
		<data>
			<parameter idref="Y.precision"/>
		</data>
	</multivariateWishartPrior>

	<!-- END Multivariate diffusion model                                        -->

	

	<!-- START Multivariate diffusion model                                      -->

	<multivariateTraitLikelihood id="X.traitLikelihood" traitName="X" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">
		<multivariateDiffusionModel idref="X.diffusionModel"/>
		<treeModel idref="treeModel"/>
		<traitParameter>
			<parameter id="leaf.X"/>
		</traitParameter>
		<conjugateRootPrior>
			<meanParameter>
				<parameter value="0.0"/>
			</meanParameter>
			<priorSampleSize>
				<parameter value="0.000001"/>
			</priorSampleSize>
		</conjugateRootPrior>
	</multivariateTraitLikelihood>
	<matrixInverse id="X.varCovar">
		<matrixParameter idref="X.precision"/>
	</matrixInverse>
	<continuousDiffusionStatistic id="X.diffusionRate">
		<multivariateTraitLikelihood idref="X.traitLikelihood"/>
	</continuousDiffusionStatistic>


	<multivariateTraitLikelihood id="Y.traitLikelihood" traitName="Y" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">
		<multivariateDiffusionModel idref="Y.diffusionModel"/>
		<treeModel idref="treeModel"/>
		<traitParameter>
			<parameter id="leaf.Y"/>
		</traitParameter>
		<conjugateRootPrior>
			<meanParameter>
				<parameter value="0.0"/>
			</meanParameter>
			<priorSampleSize>
				<parameter value="0.000001"/>
			</priorSampleSize>
		</conjugateRootPrior>
	</multivariateTraitLikelihood>
	<matrixInverse id="Y.varCovar">
		<matrixParameter idref="Y.precision"/>
	</matrixInverse>
	<continuousDiffusionStatistic id="Y.diffusionRate">
		<multivariateTraitLikelihood idref="Y.traitLikelihood"/>
	</continuousDiffusionStatistic>

	<!-- END Multivariate diffusion model                                        -->


	<!-- Define operators                                                        -->
	<operators id="operators" optimizationSchedule="log">

		<!-- START Multivariate diffusion model                                      -->
		<precisionGibbsOperator weight="1">
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			<multivariateWishartPrior idref="X.precisionPrior"/>
		</precisionGibbsOperator>
		<precisionGibbsOperator weight="1">
			<multivariateTraitLikelihood idref="Y.traitLikelihood"/>
			<multivariateWishartPrior idref="Y.precisionPrior"/>
		</precisionGibbsOperator>

		<!-- END Multivariate diffusion model                                        -->

	</operators>
	

	<!-- Define MCMC                                                             -->
	<mcmc id="mcmc" chainLength="1000000" autoOptimize="true" operatorAnalysis=""" +'"beast'+str(i)+'.ops.txt"'+""">
		<joint id="joint">
			<prior id="prior">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateWishartPrior idref="X.precisionPrior"/>
				<multivariateWishartPrior idref="Y.precisionPrior"/>

				<!-- END Multivariate diffusion model                                        -->

			</prior>
			<likelihood id="likelihood">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateTraitLikelihood idref="X.traitLikelihood"/>
				<multivariateTraitLikelihood idref="Y.traitLikelihood"/>

				<!-- END Multivariate diffusion model                                        -->

			</likelihood>
		</joint>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="1000">
			<column label="Joint" dp="4" width="12">
				<joint idref="joint"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="age(root)" sf="6" width="12">
				<tmrcaStatistic idref="age(root)"/>
			</column>
			
		</log>

		<!-- write log to file                                                       -->
		<log id="fileLog" logEvery="1000" fileName="""+'"beast'+str(i)+'.log.txt"'+""" overwrite="false">
			<joint idref="joint"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="treeModel.rootHeight"/>
			<tmrcaStatistic idref="age(root)"/>
			<treeLengthStatistic idref="treeLength"/>
			

			<!-- START Multivariate diffusion model                                      -->
			<matrixParameter idref="X.precision"/>
			<matrixInverse idref="X.varCovar"/>
			<continuousDiffusionStatistic idref="X.diffusionRate"/>
			<matrixParameter idref="Y.precision"/>
			<matrixInverse idref="Y.varCovar"/>
			<continuousDiffusionStatistic idref="Y.diffusionRate"/>

			<!-- END Multivariate diffusion model                                        -->

			<!-- START Multivariate diffusion model                                      -->
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			<multivariateTraitLikelihood idref="Y.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

			
			
		</log>

		<!-- write tree log to file                                                  -->
		<logTree id="treeFileLog" logEvery="1000" nexusFormat="true" fileName="""+'"beast'+str(i)+'.trees.txt"'""" sortTranslationTable="true">
			<treeModel idref="treeModel"/>
			
			<joint idref="joint"/>

			<!-- START Ancestral state reconstruction                                    -->
			<trait name="X" tag="X">
				<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			</trait>
			<trait name="Y" tag="Y">
				<multivariateTraitLikelihood idref="Y.traitLikelihood"/>
			</trait>

			<!-- END Ancestral state reconstruction                                      -->


			<!-- START Multivariate diffusion model                                      -->
			<multivariateDiffusionModel idref="X.diffusionModel"/>
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			<multivariateDiffusionModel idref="Y.diffusionModel"/>
			<multivariateTraitLikelihood idref="Y.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

		</logTree>
	</mcmc>
	
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>\n""")    
        
    file.write('</beast>\n')
        
    file.close()
    
def write_BEAST_xml_dim_1(t, d, i):
    file = open("output8/beast"+str(i)+".xml","w")
    file.write('<?xml version="1.0" standalone="yes"?>\n')
    file.write('<beast version="1.10.4">\n')
    file.write('\t<taxa id="taxa">\n')
    for tax in d:
        file.write('\t\t<taxon id="'+tax.label+'">\n')
        file.write('\t\t\t<date value="'+str(t.find_node_for_taxon(tax).time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).X)+'\n')
        file.write('\t\t\t</attr>\n')        
        ##perhaps not needed?
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).X)+'\n')
        file.write('\t\t\t</attr>\n')

        file.write('\t\t</taxon>\n')   
    file.write('\t</taxa>\n')  
    
    file.write('\t<newick id="startingTree">\n')
    file.write('\t\t'+t.as_string(schema="newick",suppress_rooting=True)+'\n')
    
    file.write('\t</newick>\n')
    
    file.write("""	<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>\n""")
    
    file.write("""	<!-- Statistic for sum of the branch lengths of the tree (tree length)       -->
	<treeLengthStatistic id="treeLength">
		<treeModel idref="treeModel"/>
	</treeLengthStatistic>

	<!-- Statistic for time of most recent common ancestor of tree               -->
	<tmrcaStatistic id="age(root)" absolute="true">
		<treeModel idref="treeModel"/>
	</tmrcaStatistic>

<!-- START Multivariate diffusion model                                      -->

	<multivariateDiffusionModel id="X.diffusionModel">
		<precisionMatrix>
			<matrixParameter id="X.precision">
				<parameter id="X.precision.col1" value="0.05"/>
			</matrixParameter>
		</precisionMatrix>
	</multivariateDiffusionModel>

	<multivariateWishartPrior id="X.precisionPrior" df="1">
		<scaleMatrix>
			<matrixParameter>
				<parameter value="1.0"/>
			</matrixParameter>
		</scaleMatrix>
		<data>
			<parameter idref="X.precision"/>
		</data>
	</multivariateWishartPrior>

	<!-- END Multivariate diffusion model                                        -->

	

	<!-- START Multivariate diffusion model                                      -->

	<multivariateTraitLikelihood id="X.traitLikelihood" traitName="X" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">
		<multivariateDiffusionModel idref="X.diffusionModel"/>
		<treeModel idref="treeModel"/>
		<traitParameter>
			<parameter id="leaf.X"/>
		</traitParameter>
		<conjugateRootPrior>
			<meanParameter>
				<parameter value="0.0"/>
			</meanParameter>
			<priorSampleSize>
				<parameter value="0.000001"/>
			</priorSampleSize>
		</conjugateRootPrior>
	</multivariateTraitLikelihood>
	<matrixInverse id="X.varCovar">
		<matrixParameter idref="X.precision"/>
	</matrixInverse>
	<continuousDiffusionStatistic id="X.diffusionRate">
		<multivariateTraitLikelihood idref="X.traitLikelihood"/>
	</continuousDiffusionStatistic>


	<!-- END Multivariate diffusion model                                        -->


	<!-- Define operators                                                        -->
	<operators id="operators" optimizationSchedule="log">

		<!-- START Multivariate diffusion model                                      -->
		<precisionGibbsOperator weight="1">
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			<multivariateWishartPrior idref="X.precisionPrior"/>
		</precisionGibbsOperator>

		<!-- END Multivariate diffusion model                                        -->

	</operators>
	

	<!-- Define MCMC                                                             -->
	<mcmc id="mcmc" chainLength="1000000" autoOptimize="true" operatorAnalysis=""" +'"beastfiles\beast'+str(i)+'.ops.txt"'+""">
		<joint id="joint">
			<prior id="prior">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateWishartPrior idref="X.precisionPrior"/>
				<multivariateWishartPrior idref="Y.precisionPrior"/>

				<!-- END Multivariate diffusion model                                        -->

			</prior>
			<likelihood id="likelihood">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateTraitLikelihood idref="X.traitLikelihood"/>
				<multivariateTraitLikelihood idref="Y.traitLikelihood"/>

				<!-- END Multivariate diffusion model                                        -->

			</likelihood>
		</joint>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="1000">
			<column label="Joint" dp="4" width="12">
				<joint idref="joint"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="age(root)" sf="6" width="12">
				<tmrcaStatistic idref="age(root)"/>
			</column>
			
		</log>

		<!-- write log to file                                                       -->
		<log id="fileLog" logEvery="1000" fileName="""+'"beastfiles\beast'+str(i)+'.log.txt"'+""" overwrite="false">
			<joint idref="joint"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="treeModel.rootHeight"/>
			<tmrcaStatistic idref="age(root)"/>
			<treeLengthStatistic idref="treeLength"/>
			

			<!-- START Multivariate diffusion model                                      -->
			<matrixParameter idref="X.precision"/>
			<matrixInverse idref="X.varCovar"/>
			<continuousDiffusionStatistic idref="X.diffusionRate"/>

			<!-- END Multivariate diffusion model                                        -->

			<!-- START Multivariate diffusion model                                      -->
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

			
			
		</log>
        
		<!-- write tree log to file                                                  -->
		<logTree id="treeFileLog" logEvery="1000" nexusFormat="true" fileName="""+'"beastfiles\beast'+str(i)+'.trees.txt"'""" sortTranslationTable="true">
			<treeModel idref="treeModel"/>
			
			<joint idref="joint"/>

			<!-- START Ancestral state reconstruction                                    -->
			<trait name="X" tag="X">
				<multivariateTraitLikelihood idref="X.traitLikelihood"/>
			</trait>

			<!-- END Ancestral state reconstruction                                      -->


			<!-- START Multivariate diffusion model                                      -->
			<multivariateDiffusionModel idref="X.diffusionModel"/>
			<multivariateTraitLikelihood idref="X.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

		</logTree>
	</mcmc>
	
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>\n""")    
        
    file.write('</beast>\n')
        
    file.close()
    
