package main

import (
	"log"
	"math"
	"math/rand"
	"nonlinear_regression"
	"os"
	"strings"
)

func generate_sample_data() {
	const DATA_NUM = 500

	wf, err := os.Create(sample_file)
	if err != nil {
		log.Fatal(err)
	}
	defer wf.Close()
	wf.WriteString("y\tx1\tx2\tx3\tx4\tx5\tx6\n")

	for range DATA_NUM {
		intercept := 200.0
		x1 := 8000 + rand.Float64()*5000
		x2 := 50 + rand.Float64()*6000
		x3 := 0.10 + rand.Float64()/10
		x4 := rand.Float64() * 124
		x5 := rand.Float64()*80000000 + 1000
		x6 := (rand.Float64()*9 + 1) / 8000
		noise := rand.NormFloat64() * 1000
		// x1 := rand.Float64() * 9000
		// x2 := rand.Float64() * 6000
		// x3 := rand.Float64() * 40
		// x4 := rand.Float64() * 5000
		// x5 := rand.Float64() * 7000
		// x6 := rand.Float64() * 8000
		// noise := rand.NormFloat64() * 5

		f2s := nonlinear_regression.FloatToString
		//x2 is unused
		y := intercept + 2*x1 + 3*math.Exp(72*x3) + 5*x4*x4 + 7*math.Sqrt(x5) + 8*1.0/x6 + noise
		// y := intercept + 2*x1 + 3*x3*x3 + 6*x4 + 7*x5 + 8*x6 //+ noise
		wf.WriteString(strings.Join([]string{f2s(y), f2s(x1), f2s(x2), f2s(x3), f2s(x4), f2s(x5), f2s(x6)}, "\t") + "\n")
	}
}
