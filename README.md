
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CloBalGame

<!-- badges: start -->
<!-- badges: end -->

**CloBalGame** is an R package for computing the *closest balanced game*
to a given TU game. It helps find a game as close as possible (in a
weighted Euclidean sense) to the original one, while satisfying the
conditions of balancedness. This package is compatible with other R
packages for cooperative game theory, such as *CoopGame*.

Sometimes, the value of a cooperative game $v$ is obtained through an
estimation procedure or expert judgment. Suppose the players agree that
it is beneficial for them to collaborate and form the grand coalition.
If the game has a non-empty core, any core imputation provides a
reasonable way to distribute the total gain among the players.

However, if the core is empty, this might indicate that the estimation
of the game was flawed, since an empty core suggests that forming the
grand coalition is not stable (some subcoalitions may have an incentive
to act independently).

In such cases, searching for the closest balanced game becomes a natural
approach. It allows us to adjust the game slightly in order to ensure a
non-empty core, while remaining as close as possible to the original
game in a weighted Euclidean sense.

In other words, it solves the following optimization problem:

$$
\min_{v^{\ast}, \ x^{\ast}} \sum_{S\subseteq N} \gamma_S  (v(S) - v^{\ast}(S))^2.
$$

subject to the following constraints:

- **Efficiency:** $\sum_{i \in N} x^{\ast}_i = v^{\ast}(N),$

- **Coalitional Rationality:**
  $\sum_{i \in S} x^{\ast}_i \ge v^{\ast}(S), \ \forall S \subsetneq N.$

Optionally, the following constraints can also be enforced:

- $v^{\ast}(S) \ge 0$ for all coalitions (`positive_game`).
- $v^{\ast}(N) = v_N$, to fix the value of the grand coalition (`vN`).
  By default, this value is set to match the original game,
  $v_N = v(N)$.

## Installation

You can install the development version of `CloBalGame` directly from
GitHub using:

``` r
# Install devtools if needed
# install.packages("devtools")

devtools::install_github("pgsegador/CloBalGame")
```

Alternatively, you can install the package manually:

- Download the .tar.gz source file from the GitHub releases or build it
  locally using devtools::build().

- Install it with:

``` r
# Replace with your actual file path
install.packages("CloBalGame_0.5.tar.gz", repos = NULL, type = "source")
```

After installation, load the package with:

``` r
library(CloBalGame)
```

## Example 1: Closest balanced game for a game with 4 players

Suppose we have measured the value of each coalition in a cooperative
game with 4 players, and obtained the following values:

``` r
v <- c(
   9.72, -3.74, -4.14, -9.90, -1.27,
   7.51, -7.80, 6.74, -5.06, 1.15,
  -0.03, -9.89, -2.68, 5.28, 2.14
)
```

If we want to display the coalitions as names of the game vector, we can
do the following:

``` r
v <- set_label_game(v)
print(v)
#>     1     2     3     4    12    13    14    23    24    34   123   124   134 
#>  9.72 -3.74 -4.14 -9.90 -1.27  7.51 -7.80  6.74 -5.06  1.15 -0.03 -9.89 -2.68 
#>   234  1234 
#>  5.28  2.14
```

Let’s check that the game is not balanced and, therefore, that there is
no core imputation that benefits all coalitions.

``` r
CoopGame::isBalancedGame(v)
#> [1] FALSE
```

We will now search for the closest balanced game. For small games
($n \leq 5$), it is recommended to use the quadratic method
(`method = "quad"`):

``` r
clobal <- closest_balanced_game(v,method='quad')
print(clobal)
#> $v_star
#>  [1]  2.264 -3.740 -4.140 -9.900 -1.270  7.510 -7.800  5.714 -6.086  0.124
#> [11] -0.030 -9.890 -2.680 -0.124  2.140
#> 
#> $x_star
#> [1]  2.264 -0.248  5.962 -5.838
```

The function returns the closest game `v_star` and a core imputation
`x_star` for it.

We can also verify that this game `v_star` has a unique core and that
`x_star` is the core element.

``` r
core <- CoopGame::coreVertices(clobal$v_star)
print(core)
#>       [,1]   [,2]  [,3]   [,4]
#> [1,] 2.264 -0.248 5.962 -5.838
```

The function also allows computing the closest **positive balanced
game** using the `positive_game` parameter.

``` r
clobal_pos <- closest_balanced_game(v,positive_game=TRUE, method='quad')
print(clobal_pos)
#> $v_star
#>  [1] 0.66 0.00 0.00 0.00 0.00 2.14 0.00 1.48 0.00 1.15 0.00 0.00 0.00 1.48 2.14
#> 
#> $x_star
#> [1]  6.600000e-01  0.000000e+00  1.480000e+00 -5.839773e-16
```

Moreover, suppose the value $v(N)$ of the grand coalition is approximate
and we do not want to impose that $v^{\ast}(N) = v(N)$. In that case, we
can call the function with `vN = NULL`, and the optimization will allow
$v^{\ast}(N)$ to vary freely.

``` r
clobal_free <- closest_balanced_game(v,vN=NULL, method='quad')
print(clobal_free)
#> $v_star
#>  [1]  5.334118 -3.740000 -4.140000 -9.900000 -1.270000  7.510000 -7.800000
#>  [8]  6.591176 -5.208823  1.001176 -0.030000 -9.890000 -2.680000  1.191765
#> [15]  6.525882
#> 
#> $x_star
#> [1]  5.3341176  0.1905882  6.4005882 -5.3994117
```

In the following table, we compare the values of these games:

| Coalition |     v | v_star | v_star_pos | v_star_free |
|:----------|------:|-------:|-----------:|------------:|
| 1         |  9.72 |   2.26 |       0.66 |        5.33 |
| 2         | -3.74 |  -3.74 |       0.00 |       -3.74 |
| 3         | -4.14 |  -4.14 |       0.00 |       -4.14 |
| 4         | -9.90 |  -9.90 |       0.00 |       -9.90 |
| 12        | -1.27 |  -1.27 |       0.00 |       -1.27 |
| 13        |  7.51 |   7.51 |       2.14 |        7.51 |
| 14        | -7.80 |  -7.80 |       0.00 |       -7.80 |
| 23        |  6.74 |   5.71 |       1.48 |        6.59 |
| 24        | -5.06 |  -6.09 |       0.00 |       -5.21 |
| 34        |  1.15 |   0.12 |       1.15 |        1.00 |
| 123       | -0.03 |  -0.03 |       0.00 |       -0.03 |
| 124       | -9.89 |  -9.89 |       0.00 |       -9.89 |
| 134       | -2.68 |  -2.68 |       0.00 |       -2.68 |
| 234       |  5.28 |  -0.12 |       1.48 |        1.19 |
| 1234      |  2.14 |   2.14 |       2.14 |        6.53 |

Comparison of different versions of the closest balanced game

We can observe that `v_star_pos` has all positive entries, and that the
value $v^{\ast}(N)$ associated with `v_star_free` does **not** match the
original $v(N)$, as it has been freely optimized without that
constraint.

## Example 2: Closest balanced game for a game with 8 players

Now, let’s consider a larger game with $n = 8$ players. We are going to
work with a distance metric based on the **Shapley value**:

$$\gamma_S = \frac{1}{\binom{n-2}{|S|-1}}.$$

When the number of players is large, the option `method = "quad"`
becomes inefficient. Instead, we suggest using the **iterative method**
with `method = "iter"`. This method has the downside that it may
occasionally get stuck in loops without converging to a solution. As the
number of players $n$ increases, it becomes less likely for the
algorithm to get stuck in a loop. This is related to the fact that the
probability that the core of $v^{\ast}$ is a singleton converges to one
as $n$ increases.

``` r
set.seed(123)  # Set seed for reproducibility

# Number of players in the game
n <- 8

# Simulation range for the values of the game (cube side)
l <- 10

# Decimal rounding parameter
round_par <- 2

# Compute the Shapley metric weight vector
gamma <- get_shapley_gamma(n)

# Whether to enforce the balanced game to be positive
positive_game <- FALSE

# Generate a random cooperative game with values in [-l, l]
len <- 2^n - 1
v <- 2 * l * runif(len) - l
v <- round(v, round_par)  # Round to the desired number of decimals

# Compute the closest balanced game using the iterative method
clobal <- closest_balanced_game(v, gamma=gamma, method = 'iter')
# print x_Star
print(clobal$x_star)
#> [1] -0.7993513  1.4609532 -0.2910149  2.0640987  1.6714484 -2.6808937 -1.8050627
#> [8]  0.9998223
# Extract the core of the resulting balanced game v_star
core <- CoopGame::coreVertices(clobal$v_star)
print(core)
#>            [,1]     [,2]       [,3]     [,4]     [,5]      [,6]      [,7]
#> [1,] -0.7993513 1.460953 -0.2910149 2.064099 1.671448 -2.680894 -1.805063
#>           [,8]
#> [1,] 0.9998223
```

We can see that the core allocation of $v^{\ast}$ coincides with the
only element of its core.

Although as this game has many players examining all the values of $v$
and $v^{\ast}$ may be too costly it is always good to remember that all
$v^{\ast}$ information is stored in $x^{\ast}$. Indeed,

$v^{\ast}(S) = \min \lbrace x^{\ast}(S),\ v(S) \rbrace .$

``` r
tol <- 10^(-8)
dif <- abs(pmin(allocation2Game(clobal$x_star),v)-clobal$v_star)>tol
print(any(dif)) # no differences
#> [1] FALSE
```

We can also compute the **Shapley value** for both the original game `v`
and the balanced approximation `v_star`:

``` r
shap_v <- CoopGame::shapleyValue(v)
print(shap_v)
#> [1] -0.2115952  1.3491667  0.6399048  1.5089286  1.0528333 -2.9387857 -1.2599048
#> [8]  0.4794524
shap_v_star <- CoopGame::shapleyValue(clobal$v_star)
print(shap_v_star)
#> [1] -0.2115952  1.3491667  0.6399048  1.5089286  1.0528333 -2.9387857 -1.2599048
#> [8]  0.4794524
```

In the following table we can compare the distributions associated with
the Shapley value of $v$, $v^{\ast}$ and the single core imputation of
$v^{\ast}$.

| Player |  Shapley_v | Shapley_v_star |     x_star |
|-------:|-----------:|---------------:|-----------:|
|      1 | -0.2115952 |     -0.2115952 | -0.7993513 |
|      2 |  1.3491667 |      1.3491667 |  1.4609532 |
|      3 |  0.6399048 |      0.6399048 | -0.2910149 |
|      4 |  1.5089286 |      1.5089286 |  2.0640987 |
|      5 |  1.0528333 |      1.0528333 |  1.6714484 |
|      6 | -2.9387857 |     -2.9387857 | -2.6808937 |
|      7 | -1.2599048 |     -1.2599048 | -1.8050627 |
|      8 |  0.4794524 |      0.4794524 |  0.9998223 |

Comparison of different allocations

We can see that the shapley value of $v$ and $v^*$ coincide in the
previous table.

We know that if the gammas are those associated with the Shapley value
then the partition that solves the following problem is precisely the
Shapley value:

$$
\begin{aligned}
& \min_{x \in \mathbb{R}^{n}} && \sum_{S \subseteq N} \gamma_S \left( v(S) - x(S) \right)^2, \\
& \text{s.a.} && x(N) = v(N).
\end{aligned}
$$

We can also compute the partitions that minimize other types of metrics
such as the Euclidean metric with the function **get_psi_allocation**:

``` r
# Euclidean distance, gamma=1
gamma <- rep(1,length(v))
psi_v <- get_psi_allocation(v,gamma=gamma)
print(psi_v)
#> [1] -0.011269531  1.179042969 -0.087675781 -0.008925781  0.866699219
#> [6] -0.731425781 -0.248769531 -0.337675781
clobal_euclidean <- closest_balanced_game(v, gamma=gamma, method = 'iter')
# print x_star
print(clobal_euclidean$x_star)
#> [1] -0.2496519  0.9901872 -0.2004693 -0.1230834  1.5480393 -0.7431310 -0.2083254
#> [8] -0.3935655
# psi_value of v_star
psi_v_star <- get_psi_allocation(clobal_euclidean$v_star,gamma=gamma)
print(psi_v_star)
#> [1] -0.011269531  1.179042969 -0.087675781 -0.008925781  0.866699219
#> [6] -0.731425781 -0.248769531 -0.337675781
```

In the same way we can compare these values when changing metrics:

| Player |      Psi_v | Psi_v_star |     x_star |
|-------:|-----------:|-----------:|-----------:|
|      1 | -0.0112695 | -0.0112695 | -0.2496519 |
|      2 |  1.1790430 |  1.1790430 |  0.9901872 |
|      3 | -0.0876758 | -0.0876758 | -0.2004693 |
|      4 | -0.0089258 | -0.0089258 | -0.1230834 |
|      5 |  0.8666992 |  0.8666992 |  1.5480393 |
|      6 | -0.7314258 | -0.7314258 | -0.7431310 |
|      7 | -0.2487695 | -0.2487695 | -0.2083254 |
|      8 | -0.3376758 | -0.3376758 | -0.3935655 |

Comparison of different allocations with Euclidean distance
