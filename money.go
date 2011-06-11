package money
/*

The package contains type Money... 

type Money struct {
	M	int64
}

...which usese a fixed-length guard for precision arithmetic: the 
int64 variable Guard (and its float64 and int-related variables Guardf
and Guardi.

ROunding is done on float64 to int64 by	the Rnd() function truncating 
at values less than (.5 + (1 / Guardf))	or greater than -(.5 + (1 / Guardf)) 
in the case of negative numbers. The Guard adds four decimal places 
of protection to rounding.

DP is the decimal precision, which can be changed in the DecimalPrecision()
function.  DP hold the places after the decimalplace in teh active money struct field M

The following functions are available 

Abs Returns the absolute value of Money
	(m *Money) Abs() *Money 
Add Adds two Money types
	(m *Money) Add(n *Money) *Money 
Black-Scholes (European put and call options)
	BS(s, k, t, r, v float64, putcall string) float64 
CMP Compounded Interest Rate
	CMP(fv, pv *Money, n float64) float64 
CNI Continuous Interest 
	(pv *Money) CNI(r float64, n int) *Money 
Cov Covariance
	Cov(x, y []float64) float64 
DecimalChange resets the package-wide decimal place (default is 2 decimal places)
	DecimalChange(d int) 
Div Divides one Money type from another
	(m *Money) Div(n *Money) *Money 
FVA Future Value of an Annuity
	(m *Money) FVA(r float64, n int) *Money 
Future Value of a growing Annuity 
	(pmt *Money) FVGA(r, g float64, n int) *Money 
FV Future Value (compound interest)
	(m *Money) FV(r float64, n int) *Money 
FVi Future Value Simple Interest 
	(m *Money) FVi(r float64, n int) (m *Money) 
Gett gets value of money truncating after DP (see Value() for no truncation)
	(m *Money) Gett() int64 
Get gets the float64 value of money (see Value() for int64)
	(m *Money) Get() float64 
Interest 
	(m *Money) Interest(r float64, n int) (m *Money) 
I Interest (of 1 value returned) for n periods
	I(r float64, n int) float64 
Ifl Interest (float64 arguments) 
	Ifl(r, n float64) float64 
Is Simple Interest (also see Compounded Interest Rate)
	(m *Money) Is(pv *Money) *Money 
Mean Average 
	Mean(a []float64) float64 
MP Mortgage Payment
	(m *Money) MP(r float64, n int) *Money 
Mul Multiplies two Money types
	(m *Money) Mul(n *Money) *Money 
Mulf Multiplies a Money with a float to return a money-stored type
	(m *Money) Mulf(f float64) *Money 
Neg Returns the negative value of Money
	(m *Money) Neg() *Money 
Pow is the power of Money
	(m *Money) Pow(r float64) *Money 
Present Value of a Series of cash flows (using non integer time periods) general case
	(m *Money) PV(fvs []Money, i []float64, ns []float64) *Money 
PVA Present Value of an Annuity (ordinary)
	(m *Money) PVA(r float64, n int) *Money 
Present Value of an Annuity Due
	(m *Money) PVAD(r float64, n int) *Money 
PVf Present Value of a single future Value (using non integer period)
	(m *Money) PVf(r, n float64) *Money 
PVGA Present Value of a growing Annuity
	(m *Money) PVGA(r, g float64, n int) *Money 
PVGAD Present Value of a growing Annuity due (start of period)
	(m *Money) PVGAD(r, g float64, n int) *Money 
PVi Present Value of a single future Value (using integer period)
	(m *Money) PVi(r float64, n int) *Money 
PVP Present Value (of an integer period)
	(m *Money) PVP(r float64, n, pd int) *Money 
R Regression 
	R(x, y []float64) (a, b, r float64) 
Round Rounds int64 remainder if greater than Round
	Rnd(r int64, trunc float64) int64 
SD Standard Deviation
	SD(a []float64) float64 
SDs Standard Deviation of a sample
	SDs(a []float64) float64 
Set sets the Money field M
	(m *Money) Set(x int64) *Money 
Setf sets a float 64 into a Money type for precision calculations
	(m *Money) Setf(f float64) *Money 
Sign returns the Sign of Money 1 if positive, -1 if negative
	(m *Money) Sign() int 
String for money type representation in basic monetary unit (DOLLARS CENTS)
	(m *Money) String() string 
Sub subtracts one Money type from another
	(m *Money) Sub(n *Money) *Money 
Value returns in int64 the value of Money (also see Gett, See Get() for float64)
	(m *Money) Value() int64 
*/

import (
	"fmt"
	"math"
)

type Money struct {
	M int64 // value of the integer64 Money
}

var (
	Guardi int     = 1000 // must be at least the size of StringD() max
	Guard  int64   = int64(Guardi)
	Guardf float64 = float64(Guardi)
	DP     int64   = 100         // for default of 2 decimal places => 10^2 (can be reset)
	DPf    float64 = float64(DP) // for default of 2 decimal places => 10^2 (can be reset)
	Round  = .5 + (1 / Guardf)
	Roundn = Round * -1
)

var (
	call string = "c"
	put  string = "p"
)

const (
	DBZ     = "Divide by zero"
	DTL     = "Decimal places too large"
	DLZ     = "Decimal places cannot be less than zero"
	INF     = "Calulcations results in infinity"
	INFN    = "Calulcations results in negative infinity"
	NAN     = "Not a Number"
	NOOR    = "Number out of range"
	OVFL    = "Overflow"
	UND     = "Undefined Number: non a number, or infinity"
	STRCONE = "String Conversion error"
	MAXDEC  = 18
)

// Abs Returns the absolute value of Money
func (m *Money) Abs() *Money {
	if m.M < 0 {
		m.Neg()
	}
	return m
}

// Add Adds two Money types
func (m *Money) Add(n *Money) *Money {
	r := m.M + n.M
	if (r^m.M)&(r^n.M) < 0 {
		panic(OVFL)
	}
	m.M = r
	return m
}

// Black-Scholes (European put and call options)
// C Theoretical call premium (non-dividend paying stock)
// c = sn(d1) - ke^(-rt)N(d2)
// d1 = ln(s/k) + (r + s^2/2)t
// d2 = d1- st^1/2
// v = stock value
// k = Stock strike price
// s = Spot price
// t = time to expire in years
// r = risk free rate
// v = volitilaty (sigma)
// e math.E 2.7183
// putcall = "c" for a call or "p" for a put 
func BS(s, k, t, r, v float64, putcall string) float64 {
	d1 := (math.Log(s/k) + ((r + (math.Pow(v, 2) / 2)) * t)) / (v * math.Sqrt(t))
	d2 := d1 - (v * math.Sqrt(t))
	if putcall == call {
		return math.Erf(d1)*s - (math.Erf(d2) * k * math.Pow(math.E, (-1*r*t)))
	}
	if putcall == put {
		return k*math.Pow(math.E, (-1*r*t)) - s + ((math.Erf(-1*d2) * k * math.Pow(math.E, (-1*r*t))) - (math.Erf(-1*d1) * s))
	}
	panic(NOOR)
}

// CMP Compounded Interest Rate
// i = (fv / pv) ^ (1/n) - 1
// pv = present value
// fv = future value 
// n - number of periods
// i = interest rate in percent per period
// returned as a decimal representation of the interest rate over the period
// can return NaN (Not a Number) on improbable input values (n = 0 pv = 0)
func CMP(fv, pv *Money, n float64) float64 {
	return Ifl(float64(fv.Div(pv).Get()), 1/n) - 1
}

// CNI Continuous Interest 
// fv = pv * e ^ (i * n)
// pv - principal or present value
// i - interest rate per period
// n - number of periods
// returned as a decimal representation of the interest rate over the period
func (pv *Money) CNI(r float64, n int) *Money {
	return pv.Mulf(Ifl(math.E, (r * float64(n))))
}

// Cov Covariance
// Cov(x,y) = SIGMA(XY) - (SIGMA(X) * SIGMA(Y))
// SIGMA a total of all of the elements of a
// n is the number of x,y data points
func Cov(x, y []float64) float64 {
	if len(x) == 0 {
		panic(NOOR)
	}
	if len(x) != len(y) {
		panic(NOOR)
	}
	xy := make([]float64, len(x))
	for i, _ := range x {
		xy[i] = x[i] * y[i]
	}
	xysl := xy[:]
	return Mean(xysl) - (Mean(x) * Mean(y))
}

// DecimalChange resets the package-wide decimal place (default is 2 decimal places)
func DecimalChange(d int) {
	if d < 0 {
		panic(DLZ)
	}
	if d > MAXDEC {
		panic(DTL)
	}
	var newDecimal int
	if d > 0 {
		newDecimal++
		for i := 0; i < d; i++ {
			newDecimal *= 10
		}
	}
	DPf = float64(newDecimal)
	DP = int64(newDecimal)
	return
}

// Div Divides one Money type from another
func (m *Money) Div(n *Money) *Money {
 	f := Guardf * DPf * float64(m.M) / float64(n.M) / Guardf
	i := int64(f)
	return m.Set(Rnd(i, f-float64(i)))
}

// FVA Future Value of an Annuity
// fv = pmt * ( ( 1 + n )^n - 1 ) / r 
// fv = future value 
// pmt (m) = payment per period
// r = interest rate in percent per period
// n = number of periods
// can return NaN (Not a Number) on improbable input values
func (m *Money) FVA(r float64, n int) *Money {
	return m.Mulf((I(r, n) - 1) / r)
}

// FVGA Future Value of a growing Annuity 
// fv = pmt * ( (1+i)^n - (1+r)^n ) / (i - r) 
//   when i = r, fv = pmt * n * ((1+i)^(n-1))
// pmt (m) = amount of each payment
// fv = Future Value
// r = interest rate in decimal percent
// g rate of growth in decimal percent
// n = periods
func (m *Money) FVGA(r, g float64, n int) *Money {
	if r == g {
		return m.Mulf(float64(n) * I(r, n-1))
	}
	return m.Mulf(I(r, n) - I(g, n)/(r-g))
}

// FV Future Value (compound interest)
// fv (m) = pv * ( 1 + i )^n
// pv = present value
// fv = future value (the maturity value)
// i = interest rate in percent per period
// n = number of periods
func (m *Money) FV(r float64, n int) *Money {
	return m.Mulf(I(r, n))
}

// FVi Future Value Simple Interest 
// fv = pv * (1 + ( i * n ) )
// fv - future value (maturity value)
// pv (m) - principal or present value
// i - interest rate per period
// n - number of periods
func (m *Money) FVi(r float64, n int) *Money {
	return m.Mulf(1 + r*float64(n))
}

// Gett gets value of money truncating after DP (see Value() for no truncation)
func (m *Money) Gett() int64 {
	return m.M / DP
}

// Get gets the float64 value of money (see Value() for int64)
func (m *Money) Get() float64 {
	return float64(m.M) / DPf
}

// Interest 
// I = pv * i * n
// pv (m) - principal or present value
// i - interest rate per period
// n - number of periods
// returned into a money type as an int64-decimal representation of the interest rate over the period
func (m *Money) Interest(r float64, n int) *Money {
	return m.Mulf(I(r, n))
}

// I Interest (of 1 value returned) for n periods
// i - interest rate per period
// n - number of periods
// returned as a decimal representation of the interest rate over the period
func I(r float64, n int) float64 {
	return math.Pow((1 + r), float64(n))
}

// Ifl Interest (float64 arguments) 
// r PLUS - interest rate per period - use (1 + r), E for ln etc.
// n - number of periods
// returned as a decimal representation of the interest rate over the period
func Ifl(r, n float64) float64 {
	return math.Pow(r, n)
}

// Is Simple Interest (also see Compounded Interest Rate)
// i = fv - pv
// pv = present value in Money
// m = future value (maturity value) in Money
// i = interest rate in percent per period returned as Money
// returned as a decimal representation of the interest rate over the period
func (m *Money) Is(pv *Money) *Money {
	n := new(Money)
	return m.Sub(pv).Div(n.Set(100))
}

// Mean Average 
// mean = SIGMA a / len(a)
// SIGMA a total of all of the elements of a
// len(a) = the number of values
func Mean(a []float64) float64 {
	lenA := float64(len(a))
	if lenA == 0 {
		panic(NOOR)
	}
	var sum float64
	for _, v := range a {
		sum += v
	}
	return sum / lenA
}

// MP Mortgage Payment
// pmt = loan * r * (1 + r)^n / ((1 + r)^n - 1)
// loan - loan amount
// i - note percent interest rate (not monthly rate)
// n - number of periods (ex 360 for a 30 year loan)
// returned as Money 
func (m *Money) MP(r float64, n int) *Money {
	return m.Setf(m.Get() * r * I(r/12, n) / (I(r/12, n) - 1))
}

// Mul Multiplies two Money types
func (m *Money) Mul(n *Money) *Money {
	return m.Set(m.M * n.M / DP)
}

// Mulf Multiplies a Money with a float to return a money-stored type
func (m *Money) Mulf(f float64) *Money {
	i := m.M * int64(f*Guardf*DPf)
	r := i/Guard/DP
	return m.Set(Rnd(r,float64(i) / Guardf / DPf - float64(r)))
}

// Neg Returns the negative value of Money
func (m *Money) Neg() *Money {
	if m.M != 0 {
		m.M *= -1
	}
	return m
}

// Pow is the power of Money
func (m *Money) Pow(r float64) *Money {
	return m.Setf(math.Pow(m.Get(), r))
}

// Present Value of a Series of cash flows (using non integer time periods) general case
// pv = SIGMA (n, t=0) [fv-sub(t) / ((1 + i) ^ n-sub(t))]
// pv = present value of Money
// SIGMA (n, t=0) Sum of n's beginning at time t = 0
// fv = future value in the array slice element of fvs
// fv-sub(t) Future Value at time t (subscript t)
// n-sub(t)  period at array slice position t
// i = interest rate in percent per period
// n = period in array slice ns
// fvs and ns must correspond, be the same len()
// does not work for for example fvs[j] < -100% and ns[j] = 0.5
func (m *Money) PV(fvs []Money, i []float64, ns []float64) *Money {
	m.Set(0)
	for j, _ := range fvs {
		m.Add(fvs[j].PVf(i[j], ns[j]))
	}
	return m
}
// PVA Present Value of an Annuity (ordinary)
// pv = pmt * (1 - ( 1 / ((1+i)^n)) / i)
// m = amount of each payment and solution returned
// n = periods
// r = interest rate in decimal percent
func (m *Money) PVA(r float64, n int) *Money {
	return m.Mulf((1 - (1 / I(r, n))) / r)
}

// Present Value of an Annuity Due
// pv = pmt * (1 - ( 1 / ((1+i)^n)) / i) * (i+1)
// pmt = amount of each payment
// n = periods
// r = interest rate in decimal percent
func (m *Money) PVAD(r float64, n int) *Money {
	return m.Mulf(((1 - (1 / I(r, n))) / r) * (r + 1))
}

// PVf Present Value of a single future Value (using non integer period)
// pv = fv * ( 1 / ( 1 + r )^n)
// pv = present value
// fv = future value 
// r = interest rate in percent for period
// n = period
// can panic on improbable input values if > -100% r or 1/2 period.
// called by finance.PV
func (m *Money) PVf(r, n float64) *Money {
	return m.Mulf(1 / Ifl(1+r, n))
}

// PVGA Present Value of a growing Annuity
// pv = pmt * (1 - ( (1+g) / ((1+r)^n)) / (r-g))
// pv = present value (returned as m)
// m = payment per period
// r = interest rate in percent per period
// n = number of periods
// g  = rate of growth
func (m *Money) PVGA(r, g float64, n int) *Money {
	return m.Mulf(((1 - ((1 + g) / I(r, n))) / (r - g)))
}

// PVGAD Present Value of a growing Annuity due (start of period)
// pv = pmt * (1 - ( (1+g) / ((1+r)^n)) / (r-g)) * (1+r)
// pv = present value 
// pmt = payment per period
// r = interest rate in percent per period
// n = number of periods
// g rate of growth
func (m *Money) PVGAD(r, g float64, n int) *Money {
	return m.Mulf(((1 - ((1 + g) / I(r, n))) / (r - g)) * (1 + r))
}

// PVi Present Value of a single future Value (using integer period)
// pv = fv * ( 1 / ( 1 + r )^n)
// pv = present value
// fv = future value 
// r = interest rate 
// n = period
func (m *Money) PVi(r float64, n int) *Money {
	return m.Mulf(1 / I(r, n))
}

// PVP Present Value (of an integer period)
// pv = fv * ( 1 / ( 1 + i )^n)
// pv = present value
// fv = future value 
// r = interest rate in percent per period
// n = number of periods
func (m *Money) PVP(r float64, n, pd int) *Money {
	return m.Mulf(1 / Ifl((1+(r/float64(pd))), float64(n*pd)))
}

// R Regression 
// slope(b) = (n * SIGMAXY - (SIGMA X)(SIGMA Y))) / (n * SIGMAX^2) - (SIGMAX)^2)
// Intercept(a) = (SIGMA Y - b(SIGMA X)) / n
// r-squared = (Cov(x,y) / SD(x) * SD(y))^2
// r-squared = s1 / (p1' * q1')
// s1 = n('XY) - ('X)('Y)
// p1 = (n('X2) -- ('X)2)^1/2
// q1 = (n('X2) -- ('X)2)^1/2
func R(x, y []float64) (a, b, r float64) {
	n := float64(len(x))
	var (
		sumX, sumY, sumXY, sumXsq, sumYsq float64
	)
	for _, v := range x {
		sumX += v
	}
	for _, v := range y {
		sumY += v
	}
	for i, v := range x {
		sumXY += v * y[i]
	}
	for _, v := range x {
		sumXsq += math.Pow(v, 2)
	}
	b = ((n * sumXY) - (sumX * sumY)) / ((n * sumXsq) - math.Pow(sumX, 2))
	a = (sumY - (b * sumX)) / n
	s1 := (n * sumXY) - (sumX * sumY)
	p1 := n*sumXsq - math.Pow(sumX, 2)
	p1sqrt := math.Sqrt(p1)
	for _, v := range y {
		sumYsq += math.Pow(v, 2)
	}
	q1 := n*sumYsq - math.Pow(sumY, 2)
	q1sqrt := math.Sqrt(q1)
	r = s1 / (p1sqrt * q1sqrt)
	return a, b, r
}

// Round Rounds int64 remainder if greater than Round (or lt Roundn- neg.)
// trunc = the remainder of the float64 calc
// r     = the result of the int64 cal
func Rnd(r int64, trunc float64) int64 {
	if trunc > 0 {
		if trunc > Round {
			r++
		}
	} else {
		if trunc < Roundn {
			r--
		}
	}
	return r
}

// SD Standard Deviation
// sd = sqrt(SIGMA ((a[i] - mean) ^ 2) / len(a))
// SIGMA a total of all of the elements of a
// a[i] is the ith elemant of a
// len(a) = the number of elements in the slice a
func SD(a []float64) float64 {
	var sum float64
	m := Mean(a)
	for _, v := range a {
		sum += math.Pow(v-m, 2)
	}
	return math.Sqrt(sum / float64(len(a)))
}

// SDs Standard Deviation of a sample
// sd = sqrt(SIGMA ((a[i] - mean) ^ 2) / (len(a)-1))
// SIGMA a total of all of the elements of a
// a[i] is the ith elemant of a
// len(a) = the number of elements in the slice a adjusted for sample
func SDs(a []float64) float64 {
	var sum float64
	m := Mean(a)
	for _, v := range a {
		sum += math.Pow(v-m, 2)
	}
	return math.Sqrt(sum / float64(len(a)-1))
}

// Set sets the Money field M
func (m *Money) Set(x int64) *Money {
	m.M = x
	return m
}

// Setf sets a float64 into a Money type for precision calculations
func (m *Money) Setf(f float64) *Money {
	fDPf := f * DPf
	r := int64(f*DPf)
	return m.Set(Rnd(r, fDPf-float64(r)))
}

// Sign returns the Sign of Money 1 if positive, -1 if negative
func (m *Money) Sign() int {
	if m.M < 0 {
		return -1
	}
	return 1
}

// String for money type representation in basic monetary unit (DOLLARS CENTS)
func (m *Money) String() string {
	return fmt.Sprintf("%d.%02d", m.Value()/DP, m.Abs().Value()%DP)
}

// Sub subtracts one Money type from another
func (m *Money) Sub(n *Money) *Money {
	r := m.M - n.M
	if (r^m.M)&^(r^n.M) < 0 {
		panic(OVFL)
	}
	m.M = r
	return m
}

// Value returns in int64 the value of Money (also see Gett(), See Get() for float64)
func (m *Money) Value() int64 {
	return m.M
}
