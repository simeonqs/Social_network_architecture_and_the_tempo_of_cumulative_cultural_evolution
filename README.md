# Social network architecture and the tempo of cumulative cultural evolution

The Python and R code needed to replicate results from the article:

```
Cantor M, Chimento M, Smeele SQ, He P, Papageorgiou D, Aplin LM, Farine DR. in review. 
Social network architecture and the tempo of cumulative cultural evolution. 
https://www.biorxiv.org/content/10.1101/2020.12.04.411934v1
```

##### Version 0.2 January 2021

------------------------------------------------

**Abstract**

The ability to build upon previous knowledge—cumulative cultural evolution—is a hallmark of human societies. While cumulative cultural evolution depends on the interaction between social systems, cognition and the environment, there is increasing evidence that cumulative cultural evolution is facilitated by larger and more structured societies. However, the relative importance of social network architecture as an additional factor shaping cumulative cultural evolution remains unclear. By simulating innovation and diffusion of cultural traits in populations with stereotyped social structures, we disentangle the relative contributions of network architecture from those of population size and connectivity. We demonstrate that while more structured networks, such as those found in multilevel societies, can promote the recombination of cultural traits into high-value products, they also hinder spread and make products more likely to go extinct. We find that transmission mechanisms are therefore critical in determining the outcomes of cumulative cultural evolution. Our results highlight the complex interaction between population size, structure and transmission mechanisms, with important implications for future research.


------------------------------------------------

**The folders contain:**

1. The R code to create the networks
2. The Python code to run the agent based models
3. The R code to run the agent based models
4. The R code to create the figures


**Data availability:**

The data from 5,000 runs of these models are available in this Github repository and in the following  Dryad repository. These data can be used to replicate the figures in the article.


```
Cantor M, Chimento M, Smeele SQ, He P, Papageorgiou D, Aplin LM, Farine DR. 2021. 
Data from: Social network architecture and the tempo of cumulative cultural evolution. Dryad, Dataset, https://doi.org/10.5061/dryad.3r2280gff
```


------------------------------------------------

**Methods**

We constructed two agent-based models to simulate cultural evolution on the different types of social network architectures. We first generated social networks with six different architectures—random, small-world, lattice, modular, modular lattice, and multilevel—capturing different levels and combinations of clustering and modularity. We expressed these network architectures in populations with different sizes and densities of connections (average degree), where all individuals in the network had the same degree. We then build two agent-based models to explore how network architecture affects cumulative culture evolution. Briefly, our models allow innovations of cultural products to take place along two cultural lineages, with the knowledge of new products being spread through social connections via two transmission mechanics: either one-to-many or one-to-one diffusion. Once a high level of product diversity has been reached in both lineages, agents can recombine each lineage’s products into one with a final higher payoff product (hereafter ‘recombination’). Finally, we compared the performance of agents arranged in the different network architecture in terms of time to cultural recombination (i.e. tempo), time to diffusion, and the diversity of cultural traits.

The two agent-based models differed in the mechanics of information transmission: one-to-many versus one-to-one diffusion pathways. In our first agent-based model (model 1), all agents were initialized with an inventory of three items from each of two lineages. In each simulation round (epoch), each focal agent was selected once, at random, and a partner randomly chosen from its social network connections. These agents combined one or two items from their inventory in proportion to their value into a triad of items. If this triad was a valid product, knowledge of that product was learned, spread immediately to all their network connections (one-to-many diffusion), and subsequently became available as an ingredient for making new products. Simulations finished once a recombination product (a triad that recombines specific products from both lineages) was first innovated. We ran 5,000 simulations for each of the network architecture types, sizes and densities of connections, recording time to achieve the recombination product (in epochs) and tracking the diversity of cultural innovations over time. An epoch was one simulation round in which each agent was selected once as a focal agent in random order.

Our second agent-based model (model 2) extended the first by changing the transmission mechanic and altering the set of valid combinations such that the model can run past the first innovation of either recombination product. Transmission of valid products now occurred between dyads of agents (one-to-one diffusion) prior to choosing items from their inventory, in contrast to the broadcast style of diffusion in model 1. Secondly, if a triad contained either recombination products, the final product was that recombination product. In the case where both recombination products were present in the triad, one was chosen as the final product at random. This allowed us to track the diffusion of recombination products beyond their innovation. We also ran 5,000 simulations for the same parameter space of model 1, recording time (in epochs) to cultural lineage recombination, as well as time to diffusion to the majority of the network.

------------------------------------------------

**Authors, maintainers and contact:**


1. Simeon Q. Smeele <ssmeele@ab.mpg.de>

*Cognitive & Cultural Ecology Research Group, Max Planck Institute of Animal Behavior, Radolfzell, Germany*

*Department of Human Behavior, Ecology and Culture, Max Planck Institute for Evolutionary Anthropology, Leipzig, Germany*

2. Michael C. Chimento <mchimento@ab.mpg.de>

*Cognitive & Cultural Ecology Research Group, Max Planck Institute of Animal Behavior, Radolfzell, Germany*

3. Damien R. Farine <dfarine@ab.mpg.de>

*Department of Collective Behaviour, Max Planck Institute of Animal Behavior, Konstanz, Germany*

*Department of Biology, University of Konstanz, Konstanz, Germany*

*Centre for the Advanced Study of Collective Behaviour, University of Konstanz, Konstanz, Germany*

*Department of Evolutionary Biology and Environmental Science, University of Zurich, Switzerland*

4. Mauricio Cantor <mcantor@ab.mpg.de>

*Department for the Ecology of Animal Societies, Max Planck Institute of Animal Behavior, Konstanz, Germany*

*Departamento de Ecologia e Zoologia, Universidade Federal de Santa Catarina, Florianópolis, Brazil*



------------------------------------------------

**Funding:**

This work was supported by the Max Planck Society, the Advanced Centre for Collective Behaviour, funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy (EXC 2117–422037984), the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 850859), a Max Planck Independent Group Leader Fellowship, a scholarship from the China Scholarship Council (No. 201706100183).


