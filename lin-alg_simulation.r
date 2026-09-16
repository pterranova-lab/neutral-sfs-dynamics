library(ggplot2)
library(Matrix)

simulate_sfs <- function(
  n,
  mutcount, 
  dt, 
  t_steps, 
  mutationrate = 1e-3, 
  growthrate,
  output_dir = "sfs_frames",
  timestamps = c()) {
  
  plot_steps <- rep(FALSE, t_steps)
  plot_steps[timestamps] <- TRUE

  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
           
  # --- 1. Initialize SFS ---

  # Randomly assign each mutation to a bin index from 1 to n-1
  random_indices <- sample(1:(n-1), size = mutcount, replace = TRUE)

  # Count how many mutations fell into each bin
  # We use factor() to ensure bins with 0 mutations are still included
  SFS <- as.numeric(table(factor(random_indices, levels = 1:(n-1))))
  print(paste0("Initial Mutational load: ", sum(SFS)))
  
  # --- 2. mutationratetation vector ---
  m <- rep(0, n - 1)
  m[1] <- 2 * n * mutationrate
  
  # --- 3. Drift matrix G ---
  G <- Matrix(0, nrow = n - 1, ncol = n - 1, sparse = TRUE)
  
  j <- 1:(n - 1)

  diag(G) <- -2 * j * (n - j) / n
  G[cbind(j[-length(j)], j[-1])] <- (j[-length(j)] + 1) * (n - j[-length(j)] - 1) / n
  G[cbind(j[-1], j[-length(j)])] <- (j[-1] - 1) * (n - j[-1] + 1) / n
  

  # --- 5. Time evolution ---
  for (t in 1:t_steps) {
    
    # --- 6. ggplot histogram (bar plot) ---
    if (plot_steps[t]) {
      df <- data.frame(
        VAF = 1:(n - 1) / n,
        Mutcount = SFS
      )

      p <- ggplot(df, aes(x = VAF, y = Mutcount)) +
        geom_col(fill = "steelblue", color = "steelblue") +
        labs(
          title = paste("SFS at time step", t),
          x = "Variant allele frequency",
          y = "Mutation count"
        ) +
        theme_minimal(base_size = 14)

      ggsave(
        filename = sprintf("%s/frame_%04d.png", output_dir, t),
        plot = p,
        width = 8,
        height = 6
      )
    }

    # --- 7. Euler update ---
    SFS <- as.numeric(SFS + (G %*% SFS + m) * dt)

    
    # Ensure non-negative
    SFS <- pmax(SFS, 0)
  }
  
  return(list(
    SFS_final = SFS,
    G = G,
    m = m
  ))
}

library(ggplot2)
library(Matrix)
library(expm)  # for matrix exponentials

simulate_sfs_analytic <- function(
  n,
  mutcount,
  mutationrate = 1e-3,
  timestamps = c(),
  output_dir = "sfs_frames"
) {
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # --- 1. Initialize SFS ---
  random_indices <- sample(1:(n-1), size = mutcount, replace = TRUE)
  SFS0 <- as.numeric(table(factor(random_indices, levels = 1:(n-1))))
  cat("Initial Mutational load:", sum(SFS0), "\n")
  
  # --- 2. Mutation input vector ---
  m <- rep(0, n-1)
  m[1] <- 2 * n * mutationrate  
  
  # --- 3. Drift matrix G (sparse) ---
  G <- Matrix(0, nrow = n-1, ncol = n-1, sparse = TRUE)
  j <- 1:(n-1)
  diag(G) <- -2 * j * (n - j) / n
  G[cbind(j[-length(j)], j[-1])] <- (j[-length(j)] + 1) * (n - j[-length(j)] - 1) / n
  G[cbind(j[-1], j[-length(j)])] <- (j[-1] - 1) * (n - j[-1] + 1) / n
  
  # --- 4. Analytic equilibrium ---
  SFS_eq <- as.numeric(-solve(G, m))  # ξ_eq = -G^-1 m
  
  # --- 5. Loop over requested timestamps ---
  for (t in timestamps) {
    
    # Analytic solution: ξ(t) = exp(G*t) * (SFS0 - SFS_eq) + SFS_eq
    SFS_t <- as.numeric(expm(G * t) %*% (SFS0 - SFS_eq) + SFS_eq)
    
    # Ensure non-negative (just in case)
    SFS_t <- pmax(SFS_t, 0)
    
    # --- 6. Prepare data for plotting ---
    df <- data.frame(
      VAF = 1:(n-1) / n,
      Mutcount = SFS_t
    )
    
    # --- 7. ggplot histogram ---
    p <- ggplot(df, aes(x = VAF, y = Mutcount)) +
      geom_col(fill = "steelblue", color = "steelblue") +
      labs(
        title = paste("SFS at time t =", t),
        x = "Variant allele frequency",
        y = "Mutation count"
      ) +
      theme_minimal(base_size = 14)
    
    ggsave(
      filename = sprintf("%s/frame_%04d.png", output_dir, t),
      plot = p,
      width = 8,
      height = 6
    )
  }
  
  return(list(
    SFS_final = SFS_eq,
    SFS_initial = SFS0,
    G = G,
    m = m
  ))
}

