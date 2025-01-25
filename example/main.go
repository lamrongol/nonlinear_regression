package main

import (
	"errors"
	"fmt"
	"io/fs"
	"math"
	"nonlinear_regression"
	"os"
)

const sample_file = "sample_data.tsv"

const output_file = "result.tsv"

func main() {
	if _, err := os.Stat(sample_file); errors.Is(err, fs.ErrNotExist) {
		fmt.Println("Generating sample data")
		generate_sample_data()
		fmt.Println("Generating Finished")
	}
	nonlinear_regression.NonlinearRegression(sample_file,
		nonlinear_regression.STEPWISE_AIC,
		nil,
		true,
		output_file, 0, 0, 0, 0, 0, 0.0, 0)

	//calculate for new data
	individual := nonlinear_regression.RestoreIndividualByFile(output_file)
	fmt.Println(individual.Calculate([]float64{8624, 1880, 0.1125, 18.80, 2330000, 0.0008250}))
	//You must set any value if unused parameters exist
	fmt.Println(individual.Calculate([]float64{500, math.NaN(), 0.12, -18, 30000, 0.00075}))
}
