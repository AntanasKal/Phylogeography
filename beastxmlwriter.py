# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 09:32:25 2019

@author: Antanas
"""

#writing BEAST xml file dependiing of the dimension
def write_BEAST_xml(t, d, i, dimension):
    #t is the tree
    #d is the DNA character matrix (currently not needed)
    #i is the index of the file
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
	<mcmc id="mcmc" chainLength="50000" autoOptimize="true" operatorAnalysis=""" +'"beastfiles\\beast'+str(i)+'.ops.txt"'+""">
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
		<log id="screenLog" logEvery="50">
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
		<log id="fileLog" logEvery="50" fileName="""+'"beastfiles\\beast'+str(i)+'.log.txt"'+""" overwrite="false">
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
		<logTree id="treeFileLog" logEvery="50" nexusFormat="true" fileName="""+'"beastfiles\\beast'+str(i)+'.trees.txt"'""" sortTranslationTable="true">
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
	<mcmc id="mcmc" chainLength="50000" autoOptimize="true" operatorAnalysis=""" +'"beastfiles\\beast'+str(i)+'.ops.txt"'+""">
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
		<log id="screenLog" logEvery="50">
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
		<log id="fileLog" logEvery="50" fileName="""+'"beastfiles\\beast'+str(i)+'.log.txt"'+""" overwrite="false">
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
		<logTree id="treeFileLog" logEvery="50" nexusFormat="true" fileName="""+'"beastfiles\\beast'+str(i)+'.trees.txt"'""" sortTranslationTable="true">
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