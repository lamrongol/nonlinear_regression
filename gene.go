package nonlinear_regression

import (
	"math"
	"math/rand"
)

type Gene interface {
	Dim() int //Degree of freedom
	Name() string

	Calc(x float64) float64
}

type GeneEnum struct {
	//Fixed values for each GeneEnum
	// code            int
	name            string
	dim             int
	minus_allowable bool
	//---------------------------

	scale_factor float64
}

func (g *GeneEnum) Name() string {
	return g.name
}

func (g *GeneEnum) Dim() int {
	return g.dim
}

func (g *GeneEnum) SetScaleFactor(s float64) {
	g.scale_factor = s
}

func (g *GeneEnum) String() string {
	if g.dim == 0 {
		return g.name
	} else if g.dim == 1 { //scale_factor doesn't exist
		return g.name + "\t"
	} else {
		return g.name + "\t" + FloatToString(g.scale_factor)
	}
}

func (g *GeneEnum) Calc(x float64) float64 {
	switch g.name {
	case "Unused":
		return 0.0 //math.NaN()
	case "Linear":
		return x
	case "Squared":
		return x * x
	case "Cubed":
		return x * x * x
	case "Exp":
		//TODO
		val := math.Exp(g.scale_factor * x)
		if math.IsInf(val, 1) {
			return PracticalMax
		}
		return val
	// case "ExpSquared":
	// 	return math.Exp(g.scale_factor * x * x)
	case "ExpMinus":
		return math.Exp(-g.scale_factor * x)
	// case "ExpMinusSquared":
	// 	return math.Exp(-g.scale_factor * x * x)
	case "Inverse":
		//TODO
		if x == 0.0 {
			return PracticalMax
		}
		return 1.0 / x
	case "Sqrt":
		return math.Sqrt(x)
	case "Log":
		//TODO
		if x == 0.0 {
			x = PracticalMin
		}
		return math.Log(x)
	case "Log1Plus":
		return math.Log1p(g.scale_factor * x)
	default:
		panic("GeneEnum.Calc(): name not found")
	}
}

func (g *GeneEnum) Clone() *GeneEnum {
	return &GeneEnum{name: g.name, dim: g.dim, minus_allowable: g.minus_allowable, scale_factor: g.scale_factor}
}

func new_gene_enum(name string, dim int, minus_allowable bool) *GeneEnum {
	return &GeneEnum{name: name, dim: dim, minus_allowable: minus_allowable}
}

var minus_allowable_gene_list = []*GeneEnum{
	new_gene_enum("Unused", 0, true),
	new_gene_enum("Linear", 1, true),
	new_gene_enum("Squared", 1, true),
	new_gene_enum("Cubed", 1, true),
	new_gene_enum("Exp", 2, true),
	// new_gene_enum("ExpSquared", 2, true),
	new_gene_enum("ExpMinus", 2, true),
	// new_gene_enum("ExpMinusSquared", 2, true),
	new_gene_enum("Inverse", 1, true),
}
var except_gene_list = []*GeneEnum{
	new_gene_enum("Sqrt", 1, false),
	new_gene_enum("Log", 1, false),
	new_gene_enum("Log1Plus", 2, false),
}

var gene_list = append(minus_allowable_gene_list, except_gene_list...)

var minus_allowable_gene_num = len(minus_allowable_gene_list)
var gene_kind_num = len(gene_list)

func RandomGene(minus_possible bool, scale_seed float64) *GeneEnum {
	var gene *GeneEnum
	if minus_possible {
		gene = minus_allowable_gene_list[rand.Intn(minus_allowable_gene_num)].Clone()
	} else {
		gene = gene_list[rand.Intn(gene_kind_num)].Clone()
	}
	if gene.Dim() == 2 {
		gene.SetScaleFactor(scale_seed/10.0 + rand.Float64()*9.9*scale_seed)
	}
	return gene
}

var gene_name_map = map[string]*GeneEnum{}

func GetGeneByName(name string) *GeneEnum {
	if len(gene_name_map) == 0 {
		//initialize
		for _, gene := range gene_list {
			gene_name_map[gene.Name()] = gene.Clone()
		}
	}

	return gene_name_map[name]
}
