#include <Rcpp.h>
using namespace Rcpp;

// This code contains a series of sanity checks for preprocessing.
// These checks are used to ensure that the input data is consistent
// Date: 2024-03-12
// Author: @marcelortizv
// Version: 0.0.0.9000


// [[Rcpp::export]]
void checkPartitionUniqueness(DataFrame dta, String idname, String partitionName) {
  CharacterVector idNames = dta[idname.get_cstring()];
  CharacterVector partitions = dta[partitionName.get_cstring()];
  
  // Create a map to hold the unique partition value for each id
  std::map<String, String> uniquePartitions;
  
  // Loop through the dataset
  for(int i = 0; i < idNames.size(); i++) {
    String currentID = idNames[i];
    String currentPartition = partitions[i];
    
    // Check if this id has been seen before
    if(uniquePartitions.find(currentID) == uniquePartitions.end()) {
      // If not, add it to the map
      uniquePartitions[currentID] = currentPartition;
    } else {
      // If yes, check if the partition is the same as before
      if(uniquePartitions[currentID] != currentPartition) {
        stop("The value of partition.name must be the same across all periods for each particular unit");
        //std::string errorMessage = "The value of " + std::string(partitionName.get_cstring()) + " must be the same across all periods for each particular unit";
        //stop(errorMessage.c_str());
      }
    }
  }
  // If we reach here, all partitions are unique by idname
}

// [[Rcpp::export]]
void checkTreatmentUniqueness(DataFrame dta, String idname, String treatName) {
  CharacterVector id = dta[idname];
  NumericVector cohort = dta[treatName];
  
  std::map<String, std::set<double>> cohortMap;
  
  // Building a map where each key is an idname and the value is a set of cohort.names
  for (int i = 0; i < id.size(); ++i) {
    cohortMap[id[i]].insert(cohort[i]);
  }
  
  // Check if any id has more than one unique cohort.name
  for (const auto &pair : cohortMap) {
    if (pair.second.size() > 1) {
      stop("The value of dname (treatment variable) must be the same across all periods for each particular unit");
    }
  }
  // If we reach here, all dname are unique by idname
}

// [[Rcpp::export]]
void checkWeightsUniqueness(DataFrame dta, String idname) {
  CharacterVector id = dta[idname];
  NumericVector weights = dta["weights"];
  
  std::map<String, std::set<double>> weightsMap;
  
  // Building a map where each key is an idname and the value is a set of weights
  for (int i = 0; i < id.size(); ++i) {
    weightsMap[id[i]].insert(weights[i]);
  }
  
  // Check if any id has more than one unique weight
  for (const auto &pair : weightsMap) {
    if (pair.second.size() > 1) {
      stop("The value of weights must be the same across all periods for each particular unit");
    }
  }
  // If we reach here, all weights are unique by idname
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
