import BasicOperations as bo
import Homogenization as homo

matrix, homogenization, loading, strainconditions, time =bo.ReadInputs("InputFile.txt")
homo.HomogeneousMaterialVEVPFiniteStrain(matrix,loading,strainconditions,time,"results.csv")