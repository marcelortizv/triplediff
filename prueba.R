seed = 123
num_ids = 100
time = 2
initial.year = 2019
treatment.year = 2020

sim_data = generate_test_panel(seed = seed,
                               num_ids = num_ids,
                               time = time,
                               initial.year = initial.year,
                               treatment.year = treatment.year)

ddd_test <- ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                data = sim_data, control.group = NULL, estMethod = "trad", learners = NULL,
                weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE)


summary(ddd_test)
