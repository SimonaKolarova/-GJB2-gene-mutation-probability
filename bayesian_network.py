import itertools
import sys

PROBS = {
    # Unconditional probabilities for pathogenic mutations in GJB2 gene
    "GJB2-mutations": {
        2: 0.00023,
        1: 0.02998,
        0: 0.96979
    },

    "Trait-expression": {
        # Probability of expressing trait given two pathogenic copies of GJB2 gene
        2: {
            True: 1,
            False: 0
        },

        # Probability of expressing trait given one pathogenic copy of GJB2 gene
        1: {
            True: 0.0004,
            False: 0.9996

        },

        # Probability of expressing trait given no pathogenic copies of GJB2 gene
        0: {
            True: 0.0008,
            False: 0.9992
        }
    },

    # Mutation probability
    "Mutation": 0.000000011
}

def calculate_probabilities(people):
    """Given a dictionary of `people` 
    calculates the probability of each person 
    to have none, one or two mutated copies of the GJB2 gene and
    to express the hearing imparement trait.
    """

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    return(probabilities)

def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Computes and returns a joint probability
    for combination of arguments provided.
    """
    # Initiate joint probability
    joint_probability = 1
    
    # Calculate joint gene and trait probabilities for parents
    for parent in people: 
        if people[parent]['mother'] != None and people[parent]['father'] != None:
            continue

        # Probability of parent having no copies of the gene
        if parent not in one_gene and parent not in two_genes:
            if parent in have_trait: # And expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][0]*PROBS["Trait-expression"][0][True])
            if parent not in have_trait: # And not expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][0]*PROBS["Trait-expression"][0][False])

        # Probability of parent having one copy of the gene
        if parent in one_gene: 
            if parent in have_trait: # And expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][1]*PROBS["Trait-expression"][1][True])
            if parent not in have_trait: # And not expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][1]*PROBS["Trait-expression"][1][False])
        
        # Probability of parent having two copies of the gene
        if parent in two_genes:
            if parent in have_trait: # And expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][2]*PROBS["Trait-expression"][2][True])
            if parent not in have_trait: # And not expressing the trait
                joint_probability *= (PROBS["GJB2-mutations"][2]*PROBS["Trait-expression"][2][False])

    # Calculate joint gene and trait probabilities for children
    for child in people: 
        if people[child]['mother'] == None or people[child]['father'] == None:
            continue
      
        parents = [people[child]['mother'], people[child]['father']]
        
        # Probability of child having no copies of the gene
        if child not in one_gene and child not in two_genes:
            
            # For each parent
            for parent in parents:
                # If parent has no copies of the gene
                if parent not in one_gene and parent not in two_genes:
                    joint_probability *= (1-PROBS["Mutation"]) 
                
                # If parent has one copy of the gene
                if parent in one_gene:
                    joint_probability *= ((1-PROBS["Mutation"]) + PROBS["Mutation"])/2

                # If parent has no copies of the gene
                if parent in two_genes:
                    joint_probability *= PROBS["Mutation"]

            # Probability of child expressing the trait
            if child in have_trait: 
                joint_probability *= PROBS["Trait-expression"][0][True]
            if child not in have_trait:
                joint_probability *= PROBS["Trait-expression"][0][False]

        # Probability of child having one copy of the gene
        if child in one_gene: 
            combination1_probability = 1 # Probability of mother passing on no copies and father passing on one copy
            combination2_probability = 1 # Probability of mother passing on one copy and father passing on no copies

            # For each parent
            for parent in parents:
    
                # If mother has no copies of the gene or father has two copies of the gene
                if (parent == parents[0] and parent not in one_gene and parent not in two_genes) or (parent == parents[1] and parent in two_genes):
                    combination1_probability *= (1-PROBS["Mutation"]) 
                    combination2_probability *= PROBS["Mutation"]
                
                # If either parent has one copy of the gene
                if parent in one_gene:
                    combination1_probability *= ((1-PROBS["Mutation"]) + PROBS["Mutation"])/2
                    combination2_probability *= ((1-PROBS["Mutation"]) + PROBS["Mutation"])/2
                
                # If mother has two copies of the gene or father has no copies of the gene
                if (parent == parents[0] and parent in two_genes) or (parent == parents[1]and parent not in one_gene and parent not in two_genes):
                    combination1_probability *= PROBS["Mutation"]
                    combination2_probability *= (1-PROBS["Mutation"]) 

            # Joint probability of combination 1 or combination2 occuring
            joint_probability *= (combination1_probability + combination2_probability)

            # Probability of child expressing the trait
            if child in have_trait: 
                joint_probability *= PROBS["Trait-expression"][1][True]
            if child not in have_trait:
                joint_probability *= PROBS["Trait-expression"][1][False]

        # Probability of child having two copies of the gene
        if child in two_genes:
            
            # For each parent
            for parent in parents:
    
                # If parent has no copies of the gene
                if parent not in one_gene and parent not in two_genes:
                    joint_probability *= PROBS["Mutation"]
                
                # If parent has one copy of the gene
                if parent in one_gene:
                    joint_probability *= ((1-PROBS["Mutation"]) + PROBS["Mutation"])/2

                # If parent has no copies of the gene
                if parent in two_genes:
                    joint_probability *= (1-PROBS["Mutation"])

            # Probability of child expressing the trait
            if child in have_trait: 
                joint_probability *= PROBS["Trait-expression"][2][True]
            if child not in have_trait:
                joint_probability *= PROBS["Trait-expression"][2][False]

    return joint_probability

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Updates "gene" and "trait" distribution `probabilities` 
    for each person using the joint probability `p` 
    for the combination of arguments (`one_gene`, `two_genes`, `have_trait`) provided.
    """
    for person in probabilities: 
        
        # Update probablities for number of genes 
        if person not in one_gene and person not in two_genes:
            probabilities[person]["gene"][0] += p 
        if person in one_gene:
            probabilities[person]["gene"][1] += p 
        if person in two_genes:
            probabilities[person]["gene"][2] += p 


        # Update probabilities for expressing the trait
        if person in have_trait:
            probabilities[person]["trait"][True] += p
        if person not in have_trait:
            probabilities[person]["trait"][False] += p     

def normalize(probabilities):
    """
    Updates `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1).
    """
    for person in probabilities:
        # Normalise probablities for number of genes 
        genes_sum = probabilities[person]["gene"][0] + probabilities[person]["gene"][1] + probabilities[person]["gene"][2]
        probabilities[person]["gene"][0] *= (100/genes_sum)
        probabilities[person]["gene"][1] *= (100/genes_sum)
        probabilities[person]["gene"][2] *= (100/genes_sum)

        # Normalise probabilities for expressing the trait
        traits_sum = probabilities[person]["trait"][True] + probabilities[person]["trait"][False]
        probabilities[person]["trait"][True] *= (100/traits_sum)
        probabilities[person]["trait"][False] *= (100/traits_sum)

        # Round probabilities
        probabilities[person]["gene"][0] = round(probabilities[person]["gene"][0],1)
        probabilities[person]["gene"][1] = round(probabilities[person]["gene"][1],1)
        probabilities[person]["gene"][2] = round(probabilities[person]["gene"][2],1)

        probabilities[person]["trait"][True] = round(probabilities[person]["trait"][True],1)
        probabilities[person]["trait"][False] = round(probabilities[person]["trait"][False],1)
