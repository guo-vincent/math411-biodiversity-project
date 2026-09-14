# math411-biodiversity-project

## 1. Write the definition of species evenness and richness. Find (or makeup) two simple examples, one with strong evenness (but not richness) and one with strong richness (but not evenness).

Richness: number of unique species.

Evenness: degree of distribution of species.

### Scenario #1 (strong richness, weak evenness):

Suppose we are dealing with an ecosystem housing species A to I.

Species A: 1
Species B: 2
Species C: 3
Species D: 4
Species E: 5
Species F: 6
Species G: 7
Species H: 8
Species I: 99999

This ecosystem is rich, with 9 distinct species. But one species overly dominates the distribution of individuals, making this ecosystem not even.

### Scenario #2 (strong evenness, weak richness):

Now suppose we are dealing with an ecosystem housing only 2 species: A and B.

Species A: 99
Species B: 99

This ecosystem is even. Both species A and species B have equal numbers of individuals. But there are only 2 species, making this ecosystem not rank high in richness.

## 2. Find P(2 Seal Salamanders chosen w/o replacement)
(6/11) * (5/10) = 3/11.

The probability is 3/11.

## 3. Find P(2 Imitator Salamanders chosen w/o replacement)
(5/11) * (4/10) = 2/11.

The probability is 2/11.

## 4. Find P(2 of the same species)
2/11 + 3/11 = 5/11.

The probability is 5/11.

## 5. Find P(2 chosen individuals differ in species)
1 - 5/11 = 6/11.

The probability is 6/11.

## 6. Problem 5, but there are now $n_{1}$ Imitator Salamanders and $n_{2}$ Seal Salamanders
P(2 Seal Salamanders chosen w/o replacement) = $\frac{n_{2}}{n_{1} + n_{2}}$ * $\frac{n_{2} - 1}{n_{1} + n_{2} - 1}$

P(2 Imitator Salamanders chosen w/o replacement) = $\frac{n_{1}}{n_{1} + n_{2}}$ * $\frac{n_{1} - 1}{n_{1} + n_{2} - 1}$

P(2 of the same species) = $\frac{n_{1}^{2} + n_{2}^{2} - n_{1} - n_{2}}{(n_{1} + n_{2})(n_{1} + n_{2}-1)}$ 

Then P(2 chosen individuals differ in species) = $1 - \frac{n_{1}^{2} + n_{2}^{2} - n_{1} - n_{2}}{(n_{1} + n_{2})(n_{1} + n_{2}-1)}$ 
= $\frac{2n_{1}n_{2}}{(n_{1} + n_{2})(n_{1} + n_{2}-1)}$

## 7. Problem 6, but we have a larger data set, where the $i^{th}$ species has $n_{i}$ individuals and there are $S$ species recorded.
Let $E$ = $\sum_{i=1}^{S}n_{i}$ be the total number of individuals.

The total number of ways to pick 2 individuals from $E$ individuals is $\binom{E}{2}$.

For species $i$, there are $n_{i}$ individuals, so the number of pairs both belonging to species $i$ is $\binom{n_{i}}{2}$.

There are in total $\sum_{i=1}^{S}\binom{n_{i}}{2}$ same species pairs.

Therefore, the number of possible different species pairs is $\binom{E}{2} - \sum_{i=1}^{S}\binom{n_{i}}{2}$.

When expressed as a probability, this becomes $\frac{\binom{E}{2} - \sum_{i=1}^{S}\binom{n_{i}}{2}}{\binom{E}{2}}$.