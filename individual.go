package nonlinear_regression

import (
	"bufio"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

type Individual struct {
	gene_num       int
	name_list      []string
	minus_possible []bool

	gene_list        []*GeneEnum
	coefficient_list []float64
	Intercept        float64

	AIC float64
	BIC float64
	R2  float64
}

func InitializeIndividual(parameter_num int, minus_possible_list []bool, scale_list []float64, name_list []string) Individual {
	var gene_list = []*GeneEnum{}
	for i := range parameter_num {
		gene := RandomGene(minus_possible_list[i], scale_list[i])
		gene_list = append(gene_list, gene)
	}
	return Individual{
		gene_num:         parameter_num,
		name_list:        name_list,
		minus_possible:   minus_possible_list,
		gene_list:        gene_list,
		coefficient_list: make([]float64, len(gene_list)),
	}
}

func (i *Individual) Clone() *Individual {
	clone := &Individual{}
	clone.gene_num = i.gene_num
	clone.name_list = make([]string, len(i.name_list))
	copy(clone.name_list, i.name_list)
	clone.minus_possible = make([]bool, len(i.minus_possible))
	copy(clone.minus_possible, i.minus_possible)
	clone.coefficient_list = make([]float64, len(i.coefficient_list))
	copy(clone.coefficient_list, i.coefficient_list)
	clone.Intercept = i.Intercept
	clone.AIC = i.AIC
	clone.BIC = i.BIC
	clone.R2 = i.R2
	for _, gene := range i.gene_list {
		clone.gene_list = append(clone.gene_list, gene.Clone())
	}
	return clone
}

func (mother *Individual) Combine(father *Individual, scale_list []float64) (*Individual, *Individual) {
	gene_num := mother.gene_num
	start_idx := rand.Intn(gene_num)
	end_idx := rand.Intn(gene_num)
	if end_idx == start_idx {
		end_idx = (start_idx + 1) % gene_num
	}
	child1 := mother.Clone()
	child2 := father.Clone()

	changed := false
	i := start_idx
	for i != end_idx {
		if child1.gene_list[i].Name() != child2.gene_list[i].Name() {
			changed = true
			child1.gene_list[i] = father.gene_list[i].Clone() //if (isPlus == null || isPlus(i) || GeneManager.allowMinus(father.genes(i)))
			child2.gene_list[i] = mother.gene_list[i].Clone() //if (isPlus == null || isPlus(i) || GeneManager.allowMinus(this.genes(i)))
		} else {
			if child1.gene_list[i].Dim() == 2 {
				changed = true
				s1 := child1.gene_list[i].scale_factor
				s2 := child2.gene_list[i].scale_factor

				child1.gene_list[i].SetScaleFactor((s1 + s2) / 2.0)
				if rand.Float32() > 0.5 {
					child2.gene_list[i].SetScaleFactor((s1 + s2))
				} else {
					child2.gene_list[i].SetScaleFactor(math.Min(s1, s2) / 2)
				}
			}
		}
		i = (i + 1) % gene_num
	}
	if !changed {
		//Mutation
		minus_allowable := child1.gene_list[i].minus_allowable
		child1.gene_list[i] = RandomGene(minus_allowable, scale_list[i])
		child2.gene_list[i] = RandomGene(minus_allowable, scale_list[i])
	}
	return child1, child2
}

func (i *Individual) String() string {
	s := "#Coefficient\tFunction\tScaling Factor(if exists)\n"
	s += "[Intercept]\t" + FloatToString(i.Intercept) + "\n"
	for idx, g := range i.gene_list {
		if i.name_list != nil {
			s += i.name_list[idx]
		}
		s += "\t"
		if g.Name() == "Unused" {
			s += "\tUnused\n"
		} else {
			s += FloatToString(i.coefficient_list[idx]) + "\t" + g.String() + "\n"
		}
	}
	s += "#AIC=" + FloatToString(i.AIC) + "\n"
	s += "#BIC=" + FloatToString(i.BIC) + "\n"
	s += "#R^2=" + FloatToString(i.R2) + "\n"

	return s
}

func RestoreIndividualByFile(tsv_file string) Individual {
	f, err := os.Open(tsv_file)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	//skip first line
	scanner.Scan()
	scanner.Text()

	scanner.Scan()
	l := scanner.Text()
	line := strings.Split(l, "\t")
	individual := Individual{}
	individual.Intercept, _ = strconv.ParseFloat(line[1], 64)
	individual.name_list = []string{}
	individual.coefficient_list = []float64{}
	individual.gene_list = []*GeneEnum{}
	individual.minus_possible = []bool{}

	for scanner.Scan() {
		l = scanner.Text()
		if l[0] == '#' {
			break
		}
		line := strings.Split(l, "\t")
		individual.name_list = append(individual.name_list, line[0])
		coe, _ := strconv.ParseFloat(line[1], 64)
		individual.coefficient_list = append(individual.coefficient_list, coe)
		gene := GetGeneByName(line[2])
		if gene.Dim() == 2 {
			scale_factor, _ := strconv.ParseFloat(line[3], 64)
			gene.SetScaleFactor(scale_factor)
		}
		individual.gene_list = append(individual.gene_list, gene)
		individual.minus_possible = append(individual.minus_possible, gene.minus_allowable)
	}
	line = strings.Split(l, "=")
	AIC, _ := strconv.ParseFloat(line[1], 64)
	individual.AIC = AIC
	line = strings.Split(scanner.Text(), "=")
	BIC, _ := strconv.ParseFloat(line[1], 64)
	individual.BIC = BIC
	line = strings.Split(scanner.Text(), "=")
	R2, _ := strconv.ParseFloat(line[1], 64)
	individual.R2 = R2

	individual.gene_num = len(individual.gene_list)
	return individual
}

func (i *Individual) Calculate(vals []float64) float64 {
	result := i.Intercept
	for idx := range i.gene_list {
		result += i.coefficient_list[idx] * i.gene_list[idx].Calc(vals[idx])
	}
	return result
}
