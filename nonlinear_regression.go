package nonlinear_regression

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"

	"github.com/lamrongol/regression"
)

func FloatToString(f float64) string {
	return strconv.FormatFloat(f, 'g', -1, 64)
}

type Evaluation int

const (
	STEPWISE_AIC Evaluation = iota
	STEPWISE_BIC
	R2
)

var PracticalMax = math.Log(math.MaxFloat32)
var PracticalMin = 1.0 / PracticalMax

// first line is
func NonlinearRegression(tsv_file string,
	evaluate Evaluation,
	minus_possible_list []bool,
	is_all_plus bool,
	record_file string,

	INDIVIDUAL_NUM int,
	TOP_SELECTION_NUM int,
	MUTATION_RATE float64,
	MIN_LOOP_COUNT int,
	MAX_LOOP_COUNT int,
	STOP_DIFF_RATE float64,
	MAX_DATA_NUM int) {
	if INDIVIDUAL_NUM == 0 {
		INDIVIDUAL_NUM = 500
	}
	if TOP_SELECTION_NUM == 0 {
		TOP_SELECTION_NUM = 5
	}
	if MUTATION_RATE == 0.0 {
		MUTATION_RATE = 0.10
	}
	// if MIN_LOOP_COUNT == 0 {
	// 	MIN_LOOP_COUNT = 3
	// }
	if MAX_LOOP_COUNT == 0 {
		MAX_LOOP_COUNT = 1000
	}
	if STOP_DIFF_RATE == 0 {
		STOP_DIFF_RATE = 0.000001
	}
	// if MAX_DATA_NUM == 0 {
	// 	MAX_DATA_NUM = 1000
	// }
	fmt.Println("target_file: ", tsv_file)
	if record_file == "" {
		base := filepath.Base(tsv_file)
		record_file = base[:len(base)-4] + "_result.tsv"
	}
	dependent_list := []float64{}
	parameter_list := [][]float64{}
	abs_sum_list := []float64{}

	file, err := os.Open(tsv_file)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	//first line
	scanner.Scan()
	line := strings.Split(scanner.Text(), "\t")
	if line[len(line)-1] == "" {
		line = line[:len(line)-1]
	}
	param_names := line[1:]
	gene_num := len(param_names)
	if minus_possible_list == nil {
		minus_possible_list = slices.Repeat([]bool{!is_all_plus}, gene_num)
	}
	for range gene_num {
		abs_sum_list = append(abs_sum_list, 0.0)
	}

	count := 0
	for scanner.Scan() {
		l := scanner.Text()
		line := strings.Split(l, "\t")
		dependent, _ := strconv.ParseFloat(line[0], 64)
		dependent_list = append(dependent_list, dependent)
		parameters := []float64{}
		for i, s := range line[1:] {
			param, _ := strconv.ParseFloat(s, 64)
			parameters = append(parameters, param)
			abs_sum_list[i] += math.Abs(param)
		}
		parameter_list = append(parameter_list, parameters)
		count++
		if MAX_DATA_NUM > 0 && count > MAX_DATA_NUM {
			break
		}
	}

	average_list := []float64{}
	scale_list := []float64{}
	for _, sum := range abs_sum_list {
		average := sum / float64(count)
		average_list = append(average_list, average)
		scale_list = append(scale_list, 1.0/average)
	}

	individuals := []*Individual{}
	for {
		individual := InitializeIndividual(gene_num, minus_possible_list, scale_list, param_names)
		////For test, force correct genes
		// if len(individuals) == 0 {
		// 	individual.gene_list = []*GeneEnum{GetGeneByName("Linear"), GetGeneByName("Unused"), GetGeneByName("Exp"), GetGeneByName("Squared"), GetGeneByName("Sqrt"), GetGeneByName("Inverse")}
		// 	individual.gene_list[2].SetScaleFactor(72)
		// }
		success := calc_evaluation(&individual, dependent_list, parameter_list, param_names)
		if success {
			individuals = append(individuals, &individual)
			if len(individuals) == INDIVIDUAL_NUM {
				break
			}
		}
	}

	var pre_best_evaluation float64
	switch evaluate {
	case R2:
		//For R2, bigger is better so reverse "<"
		sort.Slice(individuals, func(i, j int) bool { return individuals[i].R2 > individuals[j].R2 })
		pre_best_evaluation = individuals[0].R2
	case STEPWISE_AIC:
		sort.Slice(individuals, func(i, j int) bool { return individuals[i].AIC < individuals[j].AIC })
		pre_best_evaluation = individuals[0].AIC
	case STEPWISE_BIC:
		sort.Slice(individuals, func(i, j int) bool { return individuals[i].BIC < individuals[j].BIC })
		pre_best_evaluation = individuals[0].BIC
	}

	loop_count := 0
	for loop_count < MAX_LOOP_COUNT {
		fmt.Println("loop_count:", loop_count)
		next_individuals := []*Individual{}
		for i := range TOP_SELECTION_NUM {
			next_individuals = append(next_individuals, individuals[i])
		}

		for i := range INDIVIDUAL_NUM - TOP_SELECTION_NUM {
			idx := i + TOP_SELECTION_NUM
			if rand.Float64() < float64(TOP_SELECTION_NUM)/float64(idx) {
				next_individuals = append(next_individuals, individuals[idx])
			}
		}

		survivor_count := len(next_individuals)
		fmt.Println("survivor_count=", survivor_count)
		for len(next_individuals) < INDIVIDUAL_NUM {
			mother_idx := rand.Intn(survivor_count)
			father_idx := rand.Intn(survivor_count)
			if father_idx == mother_idx {
				father_idx = (father_idx + 1) % survivor_count
			}

			child1, child2 := next_individuals[mother_idx].Combine(next_individuals[father_idx], scale_list)
			next_individuals = append(next_individuals, child1)
			if len(next_individuals) < INDIVIDUAL_NUM {
				next_individuals = append(next_individuals, child2)
			}
		}
		//mutation
		for idx, i := range next_individuals {
			if idx < TOP_SELECTION_NUM {
				continue
			}

			if rand.Float64() < MUTATION_RATE {
				mutation_idx := rand.Intn(gene_num)
				gene := RandomGene(i.minus_possible[mutation_idx], scale_list[mutation_idx])
				i.gene_list[mutation_idx] = gene
			}
		}

		for _, i := range next_individuals {
			calc_evaluation(i, dependent_list, parameter_list, param_names)
		}
		var best_evaluation float64
		switch evaluate {
		case R2:
			//For R2, bigger is better so reverse "<"
			sort.Slice(next_individuals, func(i, j int) bool { return next_individuals[i].R2 > next_individuals[j].R2 })
			best_evaluation = next_individuals[0].R2
		case STEPWISE_AIC:
			sort.Slice(next_individuals, func(i, j int) bool { return next_individuals[i].AIC < next_individuals[j].AIC })
			best_evaluation = next_individuals[0].AIC
		case STEPWISE_BIC:
			sort.Slice(next_individuals, func(i, j int) bool { return next_individuals[i].BIC < next_individuals[j].BIC })
			best_evaluation = next_individuals[0].BIC
		}

		best := next_individuals[0]
		fmt.Println(best)

		diff_rate := math.Abs(1.0 - best_evaluation/pre_best_evaluation)
		fmt.Println("diff_rate=", diff_rate)
		switch evaluate {
		case R2:
			fmt.Println("first=", next_individuals[0].R2)
			fmt.Println("second=", next_individuals[1].R2)
		case STEPWISE_AIC:
			fmt.Println("first=", next_individuals[0].AIC)
			fmt.Println("second=", next_individuals[1].AIC)
		case STEPWISE_BIC:
			fmt.Println("first=", next_individuals[0].BIC)
			fmt.Println("second=", next_individuals[1].BIC)
		}

		individuals = next_individuals
		if loop_count > MIN_LOOP_COUNT && diff_rate < STOP_DIFF_RATE {
			// term_avg := 0.0
			// term_cnt := 0
			// avg_values := make([]float64, gene_num)
			// for i := range gene_num {
			// 	if best.gene_list[i].Name() == "Unused" {
			// 		continue
			// 	}
			// 	term := best.coefficient_list[i] * best.gene_list[i].Calc(average_list[i])
			// 	avg_values[i] = term
			// 	term_avg += term
			// 	term_cnt++
			// }
			// term_avg /= float64(term_cnt)

			// change_exist := false
			// for i := range gene_num {
			// 	if best.gene_list[i].Name() == "Unused" {
			// 		continue
			// 	}
			// 	if avg_values[i] < term_avg/float64(1000) {
			// 		alt := best.Clone()
			// 		alt.gene_list[i] = GetGeneByName("Unused")
			// 		individuals[TOP_SELECTION_NUM-1] = individuals[1]
			// 		individuals[1] = alt
			// 		change_exist = true
			// 		break
			// 	}
			// }
			// if change_exist {
			// 	continue
			// }

			break
		}
		//println(preBestEvaluation, bestEvaluation, Math.abs(preBestEvaluation - bestEvaluation), Math.abs(preBestEvaluation - bestEvaluation) / preBestEvaluation)

		pre_best_evaluation = best_evaluation
		loop_count++
		fmt.Println("------------------------------------------------------------------------")
	}

	best := individuals[0]

	wf, err := os.Create(record_file)
	if err != nil {
		log.Fatal(err)
	}
	defer wf.Close()
	wf.WriteString(best.String())

	fmt.Println("output:", record_file)
}

func calc_evaluation(individual *Individual, dependent_list []float64, parameter_list [][]float64, param_names []string) bool {
	data_num := len(dependent_list)
	r := new(regression.Regression)
	r.SetObserved("Profit")
	use_count := 0
	idx_relation := []int{}
	for i, name := range param_names {
		if individual.gene_list[i].Name() == "Unused" {
			continue
		}
		r.SetVar(i, name)
		idx_relation = append(idx_relation, i)
		use_count++
	}

	for i, parameters := range parameter_list {
		processed_parameter_list := []float64{}
		for j, param := range parameters {
			if individual.gene_list[j].Name() == "Unused" {
				continue
			}
			val := individual.gene_list[j].Calc(param)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				fmt.Println("Invalid Value:", val, "Gene:", individual.gene_list[j], "Raw Value:", param)
				panic("Cause may be minus_possible list is invalid or program itself has bug")
			}
			processed_parameter_list = append(processed_parameter_list, val)
		}
		r.Train(regression.DataPoint(dependent_list[i], processed_parameter_list))
	}
	r.Run()

	for i, c := range r.GetCoeffs() {
		//TODO
		if math.IsNaN(c) || math.IsInf(c, 0) {
			individual.R2 = 0.0
			individual.AIC = math.MaxFloat64
			individual.BIC = math.MaxFloat64
			// fmt.Println("Failed:", individual.gene_list)
			// fmt.Println("Coefficients are invalid.")
			return false
		}
		// fmt.Println("Succeeded:", individual.gene_list)

		if i == 0 {
			individual.Intercept = c
			continue
		}
		individual.coefficient_list[idx_relation[i-1]] = c
	}

	dim_sum := 0
	for _, gene := range individual.gene_list {
		dim_sum += gene.Dim()
	}

	//@see https://qiita.com/WolfMoon/items/6164c09b93ca043690b3
	d := float64(data_num) * math.Log(r.GetResidualSumOfSquares()/float64(data_num))
	individual.AIC = d + 2*float64(dim_sum+1)
	individual.BIC = d + math.Log(float64(data_num))*float64(dim_sum+1)
	individual.R2 = r.R2

	// fmt.Printf("Regression formula:\n%v\n", r.Formula)
	//fmt.Printf("Regression:\n%s\n", r)
	return true
}
