package pdc

import (
	"math/big"

	bls "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"golang.org/x/exp/rand"
)

type OrutaParty struct {
	sKey fr.Element   //用户私钥
	pKey bls.G2Affine //用户公钥
}

type OrutaParams struct {
	// Generators
	g1  bls.G1Jac
	g2  bls.G2Jac
	g1a bls.G1Affine
	g2a bls.G2Affine
}

type Oruta struct {
	n       int            // 文件数量
	k       int            // sector
	s       int            // 用户数量
	signers []OrutaParty   // 用户列表
	PP      OrutaParams    // parameters
	pak     []bls.G1Affine //public aggregate key
}

type OrutaProof struct {
	lambda []bls.G1Affine
	mu     []fr.Element
	phi    []bls.G1Affine
}

func OrutaSetup() OrutaParams {
	g1, g2, g1a, g2a := bls.Generators()
	return OrutaParams{
		g1:  g1,
		g2:  g2,
		g1a: g1a,
		g2a: g2a,
	}
}

func NewOruta(k int, s int, n int, PP OrutaParams) Oruta {
	o := Oruta{
		n:   n,
		k:   k,
		s:   s,
		PP:  PP,
		pak: make([]bls.G1Affine, k),
	}
	o.OrutaKeyGen()
	return o
}

func (o *Oruta) OrutaKeyGen() {
	parties := make([]OrutaParty, o.s)
	sKeys := make([]fr.Element, o.s)

	for i := 0; i < o.s; i++ {
		sKeys[i].SetRandom()
		gx := new(bls.G2Jac).ScalarMultiplication(&o.PP.g2, sKeys[i].BigInt(&big.Int{}))
		parties[i] = OrutaParty{
			sKey: sKeys[i],
			pKey: *new(bls.G2Affine).FromJacobian(gx),
		}
	}

	for i := 0; i < o.k; i++ {
		var r fr.Element
		r.SetRandom()
		pakJac := *new(bls.G1Jac).ScalarMultiplication(&o.PP.g1, r.BigInt(&big.Int{}))
		o.pak[i] = *new(bls.G1Affine).FromJacobian(&pakJac)
	}

	o.signers = parties
}

func (o *Oruta) OrutaSign(m [][]fr.Element) [][]bls.G1Affine {

	sigma := make([][]bls.G1Affine, o.n)
	for i := range sigma {
		sigma[i] = make([]bls.G1Affine, o.s)
	}

	for i := 0; i < o.n; i++ {
		//signer
		π := rand.Intn(o.s)

		//beta
		var beta bls.G1Jac
		Hid, _ := bls.HashToG1(intToBytes(i), []byte{})
		HidJac := *new(bls.G1Jac).FromAffine(&Hid)
		beta.AddAssign(&HidJac)

		for j := 0; j < o.k; j++ {
			etaJac := *new(bls.G1Jac).FromAffine(&o.pak[j])
			etam := *new(bls.G1Jac).ScalarMultiplication(&etaJac, m[i][j].BigInt(&big.Int{}))
			beta.AddAssign(&etam)
		}

		//sigma
		var pka bls.G2Jac
		for j := 0; j < o.s; j++ {
			if j == π {
				continue
			} else {
				var a fr.Element
				a.SetRandom()
				pkJac := *new(bls.G2Jac).FromAffine(&o.signers[j].pKey)
				pkaJac := *new(bls.G2Jac).ScalarMultiplication(&pkJac, a.BigInt(&big.Int{}))
				pka.AddAssign(&pkaJac)
				sigma_i := *new(bls.G1Jac).ScalarMultiplication(&o.PP.g1, a.BigInt(&big.Int{}))
				sigma[i][j] = *new(bls.G1Affine).FromJacobian(&sigma_i)
			}
		}
		pkaAff := *new(bls.G2Affine).FromJacobian(&pka)
		psipka := Psi(pkaAff) //未实现同构计算,仅用于模拟
		psipka.Neg(&psipka)
		psipakJac := *new(bls.G1Jac).FromAffine(&psipka)
		sigma_pi := beta
		sigma_pi.AddAssign(&psipakJac)
		var skInv fr.Element
		skInv.Inverse(&o.signers[π].sKey)
		sigma_pi.ScalarMultiplication(&sigma_pi, skInv.BigInt(&big.Int{}))
		sigma[i][π] = *new(bls.G1Affine).FromJacobian(&sigma_pi)
	}
	return sigma
}

func (o *Oruta) OrutaChalGen(c int) Chal {
	chal := Chal{
		index: make([]int, c),
		v:     make([]fr.Element, c),
	}
	chal.index = getRandomIndices(o.n, c)

	//v
	for i := 0; i < c; i++ {
		chal.v[i].SetRandom()
	}
	return chal
}

func (o *Oruta) OrutaProofGen(c int, m [][]fr.Element, tag [][]bls.G1Affine, chal Chal) OrutaProof {
	proof := OrutaProof{
		lambda: make([]bls.G1Affine, o.k),
		mu:     make([]fr.Element, o.k),
		phi:    make([]bls.G1Affine, o.s),
	}
	//lambda,mu
	for i := 0; i < o.k; i++ {
		var tau fr.Element
		tau.SetRandom()
		eta := *new(bls.G1Jac).FromAffine(&o.pak[i])
		lambdaJac := *new(bls.G1Jac).ScalarMultiplication(&eta, tau.BigInt(&big.Int{}))
		proof.lambda[i] = *new(bls.G1Affine).FromJacobian(&lambdaJac)

		var tauh fr.Element
		hlambda := hashG1ToFr(proof.lambda[i])
		tauh.Mul(&tau, &hlambda)
		proof.mu[i].Add(&proof.mu[i], &tauh)
	}
	for i := 0; i < c; i++ {
		index := chal.index[i]
		vi := chal.v[i]

		for j := 0; j < o.k; j++ {
			var ymterm fr.Element
			ymterm.Mul(&vi, &m[index][j])
			proof.mu[j].Add(&proof.mu[j], &ymterm)
		}
	}

	//phi
	sigmayJac := make([]bls.G1Jac, o.s)
	for i := 0; i < c; i++ {
		index := chal.index[i]
		vi := chal.v[i]

		for j := 0; j < o.s; j++ {
			sigmaJac := *new(bls.G1Jac).FromAffine(&tag[index][j])
			sigmayterm := *new(bls.G1Jac).ScalarMultiplication(&sigmaJac, vi.BigInt(&big.Int{}))
			sigmayJac[j].AddAssign(&sigmayterm)
		}
	}

	for i := 0; i < o.s; i++ {
		proof.phi[i] = *new(bls.G1Affine).FromJacobian(&sigmayJac[i])
	}
	return proof
}

func (o *Oruta) OrutaProofVerify(chal Chal, proof OrutaProof) bool {
	c := len(chal.index)

	//left g1
	var leftg1 bls.G1Jac
	for i := 0; i < c; i++ {
		index := chal.index[i]
		vi := chal.v[i]
		var hidyterm bls.G1Jac
		Hid, _ := bls.HashToG1(intToBytes(index), []byte{})
		HidJac := *new(bls.G1Jac).FromAffine(&Hid)
		hidyterm.ScalarMultiplication(&HidJac, vi.BigInt(&big.Int{}))
		leftg1.AddAssign(&hidyterm)
	}
	for i := 0; i < o.k; i++ {
		var etamuterm bls.G1Jac
		etaJac := *new(bls.G1Jac).FromAffine(&o.pak[i])
		etamuterm.ScalarMultiplication(&etaJac, proof.mu[i].BigInt(&big.Int{}))
		leftg1.AddAssign(&etamuterm)
	}
	leftg1Aff := *new(bls.G1Affine).FromJacobian(&leftg1)
	lhs, _ := bls.Pair([]bls.G1Affine{leftg1Aff}, []bls.G2Affine{o.PP.g2a})

	//ProdPair
	var ProdPair bls.GT
	var e3left bls.G1Jac
	for i := 0; i < o.s; i++ {
		pairterm, _ := bls.Pair([]bls.G1Affine{proof.phi[i]}, []bls.G2Affine{o.signers[i].pKey})
		ProdPair.Mul(&ProdPair, &pairterm)
	}
	for i := 0; i < o.k; i++ {
		var lambdahterm bls.G1Jac
		lambdaJac := *new(bls.G1Jac).FromAffine(&proof.lambda[i])
		hlambda := hashG1ToFr(proof.lambda[i])
		lambdahterm.ScalarMultiplication(&lambdaJac, hlambda.BigInt(&big.Int{}))
		e3left.AddAssign(&lambdahterm)
	}
	e3leftAff := *new(bls.G1Affine).FromJacobian(&e3left)
	rhs, _ := bls.Pair([]bls.G1Affine{e3leftAff}, []bls.G2Affine{o.PP.g2a})
	ProdPair.Mul(&ProdPair, &rhs)
	res := lhs.Equal(&ProdPair)
	if res {
		return true
	} else {
		return false
	}
}
