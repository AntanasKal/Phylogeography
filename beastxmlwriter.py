# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 09:32:25 2019

@author: Antanas
"""
import random

#writing BEAST xml file dependiing of the dimension
def write_BEAST_xml(t, i, dimension, mcmc, log_every, beast_input_string="output/beast_input/beast", beast_output_string="output/beast_output/beast"):
    #t is the tree
    #d is the DNA character matrix (currently not needed)
    #i is the index of the file
    if dimension==2:
        write_BEAST_xml_dim_2(t, i, mcmc, log_every, beast_input_string, beast_output_string)
    else:
        write_BEAST_xml_dim_1(t, i, mcmc, log_every, beast_input_string, beast_output_string)
        
    
def write_BEAST_xml_dim_1(t, i, mcmc, log_every, beast_input_string, beast_output_string):
    file = open(beast_input_string+str(i)+".xml","w")
    file.write('<?xml version="1.0" standalone="yes"?>\n')
    file.write('<beast version="1.10.4">\n')
    file.write('\t<taxa id="taxa">\n')
    for leaf in t.leaf_node_iter():
        file.write('\t\t<taxon id="'+leaf.taxon.label+'">\n')
        file.write('\t\t\t<date value="'+str(leaf.time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(leaf.X)+'\n')
        file.write('\t\t\t</attr>\n')        
        ##perhaps not needed?
#        file.write('\t\t\t<attr name="X">\n')
#        file.write('\t\t\t\t'+str(t.find_node_for_taxon(tax).X)+'\n')
#        file.write('\t\t\t</attr>\n')

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
	<mcmc id="mcmc" chainLength="""+'"'+str(mcmc)+'"'+""" autoOptimize="true" operatorAnalysis=""" +'"'+beast_output_string+str(i)+'.ops.txt"'+""">
		<joint id="joint">
			<prior id="prior">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateWishartPrior idref="X.precisionPrior"/>

				<!-- END Multivariate diffusion model                                        -->

			</prior>
			<likelihood id="likelihood">
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateTraitLikelihood idref="X.traitLikelihood"/>

				<!-- END Multivariate diffusion model                                        -->

			</likelihood>
		</joint>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="""+'"'+str(log_every)+'"'+""">
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
		<log id="fileLog" logEvery="""+'"'+str(log_every)+'"'+""" fileName="""+'"'+beast_output_string+str(i)+'.log.txt"'+""" overwrite="false">
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
		<logTree id="treeFileLog" logEvery="""+'"'+str(log_every)+'"'+""" nexusFormat="true" fileName="""+'"'+beast_output_string+str(i)+'.trees.txt"'""" sortTranslationTable="true">
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
    
    
def write_BEAST_xml_dim_2(t, i, mcmc, log_every, beast_input_string, beast_output_string):
    file = open(beast_input_string+str(i)+".xml","w")
    file.write('<?xml version="1.0" standalone="yes"?>\n')
    file.write('<beast version="1.10.4">\n')
    file.write('\t<taxa id="taxa">\n')
    for leaf in t.leaf_node_iter():
        file.write('\t\t<taxon id="'+leaf.taxon.label+'">\n')
        file.write('\t\t\t<date value="'+str(leaf.time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(leaf.X)+'\n')
        file.write('\t\t\t</attr>\n')
        file.write('\t\t\t<attr name="Y">\n')
        file.write('\t\t\t\t'+str(leaf.Y)+'\n')
        
        
        file.write('\t\t\t</attr>\n')
        
        ##perhaps not needed?
        file.write('\t\t\t<attr name="location">\n')
        file.write('\t\t\t\t'+str(leaf.X)+" " +str(leaf.Y)+'\n')
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

	<multivariateDiffusionModel id="location.diffusionModel">
		<precisionMatrix>
			<matrixParameter id="location.precision">
				<parameter id="location.precision.col1" value="0.05 0.002"/>
				<parameter id="location.precision.col2" value="0.002 0.05"/>
			</matrixParameter>
		</precisionMatrix>
	</multivariateDiffusionModel>

	<multivariateWishartPrior id="location.precisionPrior" df="2">
		<scaleMatrix>
			<matrixParameter>
				<parameter value="1.0 0.0"/>
				<parameter value="0.0 1.0"/>
			</matrixParameter>
		</scaleMatrix>
		<data>
			<parameter idref="location.precision"/>
		</data>
	</multivariateWishartPrior>

	<!-- END Multivariate diffusion model                                        -->

	

	<!-- START Multivariate diffusion model                                      -->

	<multivariateTraitLikelihood id="location.traitLikelihood" traitName="location" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">
		<multivariateDiffusionModel idref="location.diffusionModel"/>
		<treeModel idref="treeModel"/>
		<traitParameter>
			<parameter id="leaf.location"/>
		</traitParameter>
		<conjugateRootPrior>
			<meanParameter>
				<parameter value="0.0 0.0"/>
			</meanParameter>
			<priorSampleSize>
				<parameter value="0.000001"/>
			</priorSampleSize>
		</conjugateRootPrior>
	</multivariateTraitLikelihood>
	<correlation id="location.correlation" dimension1="1" dimension2="2">
		<matrixParameter idref="location.precision"/>
	</correlation>
	<matrixInverse id="location.varCovar">
		<matrixParameter idref="location.precision"/>
	</matrixInverse>
	<continuousDiffusionStatistic id="location.diffusionRate">
		<multivariateTraitLikelihood idref="location.traitLikelihood"/>
	</continuousDiffusionStatistic>

	<!-- END Multivariate diffusion model                                        -->

	<!-- Define operators                                                        -->
	<operators id="operators" optimizationSchedule="log">

		<!-- START Multivariate diffusion model                                      -->
		<precisionGibbsOperator weight="2">
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>
			<multivariateWishartPrior idref="location.precisionPrior"/>
		</precisionGibbsOperator>

		<!-- END Multivariate diffusion model                                        -->

	</operators>
	

	<!-- Define MCMC                                                             -->
	<mcmc id="mcmc" chainLength="""+'"'+str(mcmc)+'"'+""" autoOptimize="true" operatorAnalysis=""" +'"'+beast_output_string+str(i)+'.ops.txt"'+""">
		<joint id="joint">
			<prior id="prior">
				
				

				<!-- START Multivariate diffusion model                                      -->
				<multivariateWishartPrior idref="location.precisionPrior"/>

				<!-- END Multivariate diffusion model                                        -->

			</prior>
			<likelihood id="likelihood">

				<!-- START Multivariate diffusion model                                      -->
				<multivariateTraitLikelihood idref="location.traitLikelihood"/>

				<!-- END Multivariate diffusion model                                        -->

			</likelihood>
		</joint>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="""+'"'+str(log_every)+'"'+""">
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
		<log id="fileLog" logEvery="""+'"'+str(log_every)+'"'+""" fileName="""+'"'+beast_output_string+str(i)+'.log.txt"'+""" overwrite="false">
			<joint idref="joint"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="treeModel.rootHeight"/>
			<tmrcaStatistic idref="age(root)"/>
			<treeLengthStatistic idref="treeLength"/>
			

			<!-- START Multivariate diffusion model                                      -->
			<matrixParameter idref="location.precision"/>
			<correlation idref="location.correlation"/>
			<matrixInverse idref="location.varCovar"/>
			<continuousDiffusionStatistic idref="location.diffusionRate"/>

			<!-- END Multivariate diffusion model                                        -->


			<!-- START Multivariate diffusion model                                      -->
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

			
			
		</log>

		<!-- write tree log to file                                                  -->
		<logTree id="treeFileLog" logEvery="""+'"'+str(log_every)+'"'+""" nexusFormat="true" fileName="""+'"'+beast_output_string+str(i)+'.trees.txt"'""" sortTranslationTable="true">
			<treeModel idref="treeModel"/>
			<joint idref="joint"/>

			<!-- START Ancestral state reconstruction                                    -->
			<trait name="location" tag="location">
				<multivariateTraitLikelihood idref="location.traitLikelihood"/>
			</trait>

			<!-- END Ancestral state reconstruction                                      -->


			<!-- START Multivariate diffusion model                                      -->
			<multivariateDiffusionModel idref="location.diffusionModel"/>
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>

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
    
    
    
def write_BEAST_xml_dim_2_old(t, d, i):
    file = open("beast_input/beast"+str(i)+".xml","w")
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
	<mcmc id="mcmc" chainLength="10000" autoOptimize="true" operatorAnalysis=""" +'"'+'beast_output\\beast'+str(i)+'.ops.txt"'+""">
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
		<log id="screenLog" logEvery="10">
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
		<log id="fileLog" logEvery="10" fileName="""+'"'+'beast_output\\beast'+str(i)+'.log.txt"'+""" overwrite="false">
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
		<logTree id="treeFileLog" logEvery="10" nexusFormat="true" fileName="""+'"'+'s\\beast'+str(i)+'.trees.txt"'""" sortTranslationTable="true">
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
    
    
    
    
    
def write_BEAST_xml_corrected(tree, sampled_t, d, i, mcmc, log_every, beast_input_string, beast_output_string, other_sample_size, seq_len):
    
    label_sample_space = []

    for leaf in tree.leaf_node_iter():
    #print(leaf.taxon.label)
        if len(d.__getitem__((leaf.taxon.label)))<seq_len:
            label_sample_space.append(leaf.taxon.label)
    #print("possible samples")
    #print(label_sample_space)
    #print("full taxon namespace")
    #print(tree.taxon_namespace)
    
    sample_set = random.sample(range(len(label_sample_space)), other_sample_size)
    
    
    
    
    
    file = open(beast_input_string+str(i)+".xml","w")
    file.write('<?xml version="1.0" standalone="yes"?>\n')
    file.write('<beast version="1.10.4">\n')
    file.write('\t<taxa id="taxa">\n')
    for leaf in sampled_t.leaf_node_iter():
        file.write('\t\t<taxon id="'+leaf.taxon.label+'">\n')
        file.write('\t\t\t<date value="'+str(leaf.time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(leaf.X)+'\n')
        file.write('\t\t\t</attr>\n')
        file.write('\t\t\t<attr name="Y">\n')
        file.write('\t\t\t\t'+str(leaf.Y)+'\n')       
        file.write('\t\t\t</attr>\n')
        
        ##perhaps not needed?
        file.write('\t\t\t<attr name="location">\n')
        file.write('\t\t\t\t'+str(leaf.X)+" " +str(leaf.Y)+'\n')
        file.write('\t\t\t</attr>\n')

        file.write('\t\t</taxon>\n')   
        
    for index in range(other_sample_size):  
        print(sample_set[index])
        print(label_sample_space[sample_set[index]])
        print(tree.taxon_namespace.has_taxon_label(label_sample_space[sample_set[index]]))
        leaf = tree.find_node_with_taxon_label(label_sample_space[sample_set[index]])
        print(leaf)
        
        file.write('\t\t<taxon id="'+leaf.taxon.label+'">\n')
        file.write('\t\t\t<date value="'+str(leaf.time)+'" direction="forwards" units="years"/>\n')
        file.write('\t\t\t<attr name="X">\n')
        file.write('\t\t\t\t'+str(leaf.X)+'\n')
        file.write('\t\t\t</attr>\n')
        file.write('\t\t\t<attr name="Y">\n')
        file.write('\t\t\t\t'+str(leaf.Y)+'\n') 
        file.write('\t\t\t</attr>\n')
        
        ##perhaps not needed?
        file.write('\t\t\t<attr name="location">\n')
        file.write('\t\t\t\t'+str(leaf.X)+" " +str(leaf.Y)+'\n')
        file.write('\t\t\t</attr>\n')

        file.write('\t\t</taxon>\n')   
    file.write('\t</taxa>\n') 
    
    file.write('\t<alignment id="alignment" dataType="nucleotide">\n')
    for leaf in sampled_t.leaf_node_iter():
        file.write("\t\t<sequence>\n")
        file.write('\t\t\t<taxon idref="'+leaf.taxon.label+'"/>\n')
        file.write('\t\t\t'+str(d.__getitem__(leaf.taxon.label))+'\n')
        file.write("\t\t</sequence>\n")  
    for index in range(other_sample_size):
        file.write("\t\t<sequence>\n")
        file.write('\t\t\t<taxon idref="'+label_sample_space[sample_set[index]]+'"/>\n') 
        file.write('\t\t\t')
        for nucleotide in range(seq_len):
            file.write('-')
        file.write('\n')
        file.write("\t\t</sequence>\n")  
    file.write('\t</alignment>\n')
    file.write("""	<patterns id="patterns" from="1" strip="false">
		<alignment idref="alignment"/>
	</patterns>
	

	<!-- A prior assumption that the population size has remained constant       -->
	<!-- throughout the time spanned by the genealogy.                           -->
	<constantSize id="constant" units="substitutions">
		<populationSize>
			<parameter id="constant.popSize" value="1.0" lower="0.0"/>
		</populationSize>
	</constantSize>
	

	<!-- Generate a random starting tree under the coalescent process            -->
	<coalescentSimulator id="startingTree">
		<taxa idref="taxa"/>
		<constantSize idref="constant"/>
	</coalescentSimulator>
	

	<!-- Generate a tree model                                                   -->
	<treeModel id="treeModel">
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
	</treeModel>

	<!-- Statistic for sum of the branch lengths of the tree (tree length)       -->
	<treeLengthStatistic id="treeLength">
		<treeModel idref="treeModel"/>
	</treeLengthStatistic>
	

	<!-- Generate a coalescent likelihood                                        -->
	<coalescentLikelihood id="coalescent">
		<model>
			<constantSize idref="constant"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>
	

	<!-- The strict clock (Uniform rates across branches)                        -->
	<strictClockBranchRates id="branchRates">
		<rate>
			<parameter id="clock.rate" value="1.0"/>
		</rate>
	</strictClockBranchRates>
	
	<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<strictClockBranchRates idref="branchRates"/>
	</rateStatistic>
	

	<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)             -->
	<HKYModel id="hky">
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<kappa>
			<parameter id="kappa" value="2.0" lower="0.0"/>
		</kappa>
	</HKYModel>

	<!-- site model                                                              -->
	<siteModel id="siteModel">
		<substitutionModel>
			<HKYModel idref="hky"/>
		</substitutionModel>
	</siteModel>

	<!--                                                                         -->
	<statistic id="mu" name="mu">
		<siteModel idref="siteModel"/>
	</statistic>
	
	

	<!-- START Multivariate diffusion model                                      -->

	<multivariateDiffusionModel id="location.diffusionModel">
		<precisionMatrix>
			<matrixParameter id="location.precision">
				<parameter id="location.precision.col1" value="0.05 0.002"/>
				<parameter id="location.precision.col2" value="0.002 0.05"/>
			</matrixParameter>
		</precisionMatrix>
	</multivariateDiffusionModel>

	<multivariateWishartPrior id="location.precisionPrior" df="2">
		<scaleMatrix>
			<matrixParameter>
				<parameter value="1.0 0.0"/>
				<parameter value="0.0 1.0"/>
			</matrixParameter>
		</scaleMatrix>
		<data>
			<parameter idref="location.precision"/>
		</data>
	</multivariateWishartPrior>

	<!-- END Multivariate diffusion model                                        -->


	<!-- Likelihood for tree given sequence data                                 -->
	<treeDataLikelihood id="treeLikelihood" useAmbiguities="false">
		<partition>
			<patterns idref="patterns"/>
			<siteModel idref="siteModel"/>
		</partition>
		<treeModel idref="treeModel"/>
		<strictClockBranchRates idref="branchRates"/>
	</treeDataLikelihood>
	

	<!-- START Multivariate diffusion model                                      -->

	<multivariateTraitLikelihood id="location.traitLikelihood" traitName="location" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">
		<multivariateDiffusionModel idref="location.diffusionModel"/>
		<treeModel idref="treeModel"/>
		<traitParameter>
			<parameter id="leaf.location"/>
		</traitParameter>
		<conjugateRootPrior>
			<meanParameter>
				<parameter value="0.0 0.0"/>
			</meanParameter>
			<priorSampleSize>
				<parameter value="0.000001"/>
			</priorSampleSize>
		</conjugateRootPrior>
	</multivariateTraitLikelihood>
	<correlation id="location.correlation" dimension1="1" dimension2="2">
		<matrixParameter idref="location.precision"/>
	</correlation>
	<matrixInverse id="location.varCovar">
		<matrixParameter idref="location.precision"/>
	</matrixInverse>
	<continuousDiffusionStatistic id="location.diffusionRate">
		<multivariateTraitLikelihood idref="location.traitLikelihood"/>
	</continuousDiffusionStatistic>

	<!-- END Multivariate diffusion model                                        -->


	<!-- Define operators                                                        -->
	<operators id="operators" optimizationSchedule="log">
		<scaleOperator scaleFactor="0.75" weight="1">
			<parameter idref="kappa"/>
		</scaleOperator>
		<deltaExchange delta="0.01" weight="1">
			<parameter idref="frequencies"/>
		</deltaExchange>
		<scaleOperator scaleFactor="0.75" scaleAll="true" ignoreBounds="true" weight="3">
			<parameter idref="treeModel.allInternalNodeHeights"/>
		</scaleOperator>
		<subtreeSlide size="1.0" gaussian="true" weight="30">
			<treeModel idref="treeModel"/>
		</subtreeSlide>
		<narrowExchange weight="30">
			<treeModel idref="treeModel"/>
		</narrowExchange>
		<wideExchange weight="3">
			<treeModel idref="treeModel"/>
		</wideExchange>
		<wilsonBalding weight="3">
			<treeModel idref="treeModel"/>
		</wilsonBalding>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="treeModel.rootHeight"/>
		</scaleOperator>
		<uniformOperator weight="30">
			<parameter idref="treeModel.internalNodeHeights"/>
		</uniformOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>

		<!-- START Multivariate diffusion model                                      -->
		<precisionGibbsOperator weight="2">
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>
			<multivariateWishartPrior idref="location.precisionPrior"/>
		</precisionGibbsOperator>

		<!-- END Multivariate diffusion model                                        -->

	</operators>
""")
   # <mcmc id="mcmc" chainLength="""+'"'+str(mcmc)+'"'+""" autoOptimize="true" operatorAnalysis=""" +'"'+beast_output_string+str(i)+'.ops.txt"'+""">
    file.write('\t<mcmc id="mcmc" chainLength="'+str(mcmc)+'" autoOptimize="true" operatorAnalysis="'+beast_output_string+str(i)+'.ops.txt'+'">\n')
    file.write("""		<joint id="joint">
			<prior id="prior">
				<logNormalPrior mu="1.0" sigma="1.25" offset="0.0">
					<parameter idref="kappa"/>
				</logNormalPrior>
				<dirichletPrior alpha="1.0" sumsTo="1.0">
					<parameter idref="frequencies"/>
				</dirichletPrior>
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
				<coalescentLikelihood idref="coalescent"/>
				
				
				<strictClockBranchRates idref="branchRates"/>

				<!-- START Multivariate diffusion model                                      -->
				<multivariateWishartPrior idref="location.precisionPrior"/>

				<!-- END Multivariate diffusion model                                        -->

			</prior>
			<likelihood id="likelihood">
				<treeDataLikelihood idref="treeLikelihood"/>

				<!-- START Multivariate diffusion model                                      -->
				<multivariateTraitLikelihood idref="location.traitLikelihood"/>

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
			<column label="rootHeight" sf="6" width="12">
				<parameter idref="treeModel.rootHeight"/>
			</column>
		</log>
""")
#		<log id="fileLog" logEvery="""+'"'+str(log_every)+'"'+""" fileName="""+'"'+beast_output_string+str(i)+'.log.txt"'+""" overwrite="false">
    file.write('\t\t<log id="fileLog" logEvery="'+str(log_every)+'" fileName="'+beast_output_string+str(i)+'log.txt" overwrite="false">\n')
    file.write("""			<joint idref="joint"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="treeModel.rootHeight"/>
			<treeLengthStatistic idref="treeLength"/>
			<parameter idref="constant.popSize"/>
			<parameter idref="kappa"/>
			<parameter idref="frequencies"/>
			<parameter idref="clock.rate"/>
			<rateStatistic idref="meanRate"/>

			<!-- START Multivariate diffusion model                                      -->
			<matrixParameter idref="location.precision"/>
			<correlation idref="location.correlation"/>
			<matrixInverse idref="location.varCovar"/>
			<continuousDiffusionStatistic idref="location.diffusionRate"/>

			<!-- END Multivariate diffusion model                                        -->

			<treeDataLikelihood idref="treeLikelihood"/>
			<strictClockBranchRates idref="branchRates"/>

			<!-- START Multivariate diffusion model                                      -->
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

			<coalescentLikelihood idref="coalescent"/>
			
		</log>

""")
    
#		<logTree id="treeFileLog" logEvery="""+'"'+str(log_every)+'"'+""" nexusFormat="true" fileName="""+'"'+beast_output_string+str(i)+'.trees.txt"'""" sortTranslationTable="true">
    file.write('\t\t<logTree id="treeFileLog" logEvery="'+str(log_every)+'" nexusFormat="true" fileName="'+beast_output_string+str(i)+'.trees.txt'+'" sortTranslationTable="true">\n')
    file.write("""			<treeModel idref="treeModel"/>
			<trait name="rate" tag="rate">
				<strictClockBranchRates idref="branchRates"/>
			</trait>
			<joint idref="joint"/>

			<!-- START Ancestral state reconstruction                                    -->
			<trait name="location" tag="location">
				<multivariateTraitLikelihood idref="location.traitLikelihood"/>
			</trait>

			<!-- END Ancestral state reconstruction                                      -->


			<!-- START Multivariate diffusion model                                      -->
			<multivariateDiffusionModel idref="location.diffusionModel"/>
			<multivariateTraitLikelihood idref="location.traitLikelihood"/>

			<!-- END Multivariate diffusion model                                        -->

		</logTree>
	</mcmc>
	
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>
	
</beast>""")
    
    
    