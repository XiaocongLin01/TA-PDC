package pdc

import (
	"crypto/sha256"
	"encoding/binary"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	bls "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr/fft"
	"golang.org/x/exp/rand"
)

type Party struct {
	sKey1   fr.Element   //用户私钥1
	sKey2   fr.Element   //用户私钥2
	pKeyAff bls.G1Affine //用户公钥
}

type IPAProof struct {
	qTau bls.G1Affine //内积论证证据1
	rTau bls.G1Affine //内积论证证据2
	pTau bls.G1Affine // 内积论证证据3
}

type Profit struct {
	cTau bls.G1Affine // 贡献向量承诺
	P    int          // 权益
	pi   IPAProof     // 内积论证证据
}

type BatchProfit struct {
	cTau []bls.G1Affine // 贡献向量承诺
	P    []int          // 权益
	pi   IPAProof       // 内积论证证据
}

type Tag struct {
	t     []bls.G1Affine //tag-1
	sigma []bls.G1Affine //tag-2
	K     []bls.G1Affine //tag-3
	xba   []fr.Element   //tag-4
	yba   []fr.Element   //tag-5
}

type Proof struct {
	Psigma bls.G1Affine
	X      fr.Element
	Y      fr.Element
	mu     []fr.Element
	theta  bls.G1Affine
}

type Chal struct {
	index []int        //挑战的块
	v     []fr.Element //随机值1
}

type Params struct {
	// Generators
	g1       bls.G1Jac
	g2       bls.G2Jac
	g1a      bls.G1Affine
	g2a      bls.G2Affine
	g1InvAff bls.G1Affine //用于加速配对验证
	g1B      bls.G1Jac
	g1Ba     bls.G1Affine
	g2Ba     bls.G2Affine
	h1       bls.G1Jac
	h1a      bls.G1Affine
	h2a      bls.G2Affine
	hTauHAff bls.G2Affine //h^tau
	// CRS
	tau       fr.Element // FIXME: To remove, only added for testing purposes
	domain    *fft.Domain
	H         []fr.Element
	g2Tau     bls.G2Affine   //g^tau
	PoT       []bls.G1Affine // [g^{tau^i}]
	PoTH      []bls.G1Affine // [h^{tau^i}
	lagHTaus  []bls.G1Affine // [Lag_i(tau)]
	lagHTausH []bls.G1Affine // [h^Lag_i(tau)]
	lag2HTaus []bls.G2Affine // [g2^Lag_i(tau)]
}

type VK struct {
	w1Tau  bls.G1Affine // 类型价值向量承诺
	w2Tau  bls.G2Affine // 类型价值向量承诺
	vHTau  bls.G2Affine //消除函数承诺
	g2NInv bls.G2Affine // g^{1/d}
	h2NInv bls.G2Affine // h^{1/d}
}

type WTS struct {
	weights []int   // 类型价值向量
	n       int     // 文件数量
	d       int     // 数据类型数量
	s       int     // 用户数量
	signers []Party // 用户列表
	PP      Params  // parameters
	vk      VK      // 公开验证密钥
}

func Setup(d int) Params {
	//生成生成元和辅助计算元素
	g1, g2, g1a, g2a := bls.Generators()

	var tau, beta, hF fr.Element
	tau.SetRandom()  // CRS陷门
	beta.SetRandom() // 多项式系数分离用的随机因子
	hF.SetRandom()   // 用于把g生成元转换为h生成元
	tauH := *new(fr.Element).Mul(&tau, &hF)

	g1B := new(bls.G1Jac).ScalarMultiplication(&g1, beta.BigInt(&big.Int{}))
	g2Ba := new(bls.G2Affine).ScalarMultiplication(&g2a, beta.BigInt(&big.Int{}))
	h1a := new(bls.G1Affine).ScalarMultiplication(&g1a, hF.BigInt(&big.Int{}))
	h2a := new(bls.G2Affine).ScalarMultiplication(&g2a, hF.BigInt(&big.Int{}))
	h1 := new(bls.G1Jac).ScalarMultiplication(&g1, hF.BigInt(&big.Int{}))

	//g2^tau
	g2Tau := new(bls.G2Jac).ScalarMultiplication(&g2, tau.BigInt(&big.Int{}))
	// h2^tau
	hTauHAff := new(bls.G2Affine).ScalarMultiplication(&g2a, tauH.BigInt(&big.Int{}))

	//选择H
	domain := GetDomain(uint64(d))
	omH := domain.Generator
	H := make([]fr.Element, d)
	H[0].SetOne()
	for i := 1; i < d; i++ {
		H[i].Mul(&omH, &H[i-1])
	}

	// 计算PoT,PoTH
	poT := make([]fr.Element, d)
	poT[0].SetOne()
	for i := 1; i < d; i++ {
		poT[i].Mul(&poT[i-1], &tau)
	}
	PoT := bls.BatchScalarMultiplicationG1(&g1a, poT)
	PoTH := bls.BatchScalarMultiplicationG1(h1a, poT)

	// 加速计算承诺
	lagH := GetAllLagAtWithOmegas(H, tau)
	lagHTaus := bls.BatchScalarMultiplicationG1(&g1a, lagH)
	lagHTausH := bls.BatchScalarMultiplicationG1(h1a, lagH)
	lag2HTaus := bls.BatchScalarMultiplicationG2(&g2a, lagH)

	return Params{
		g1:        g1,
		g2:        g2,
		g1a:       g1a, //生成元
		g2a:       g2a, //生成元
		g1B:       *g1B,
		g1Ba:      *new(bls.G1Affine).FromJacobian(g1B),
		g2Ba:      *g2Ba,
		g1InvAff:  *new(bls.G1Affine).FromJacobian(new(bls.G1Jac).Neg(&g1)),
		h1:        *h1,
		h1a:       *h1a, //生成元
		h2a:       *h2a, //生成元
		hTauHAff:  *hTauHAff,
		domain:    domain,
		H:         H,
		tau:       tau,
		g2Tau:     *new(bls.G2Affine).FromJacobian(g2Tau),
		PoT:       PoT,
		PoTH:      PoTH,
		lagHTaus:  lagHTaus,
		lagHTausH: lagHTausH,
		lag2HTaus: lag2HTaus,
	}
}

func NewWTS(d int, s int, n int, weights []int, PP Params) WTS {
	w := WTS{
		d:       d,
		s:       s,
		n:       n,
		weights: weights,
		PP:      PP,
	}
	w.keyGen()
	return w
}

func (w *WTS) keyGen() {
	// 生成用户公私钥
	parties := make([]Party, w.s)
	sKeys1 := make([]fr.Element, w.s)
	sKeys2 := make([]fr.Element, w.s)

	for i := 0; i < w.s; i++ {
		sKeys1[i].SetRandom()
		sKeys2[i].SetRandom()
		gx := new(bls.G1Jac).ScalarMultiplication(&w.PP.g1, sKeys1[i].BigInt(&big.Int{}))
		hy := new(bls.G1Jac).ScalarMultiplication(&w.PP.h1, sKeys2[i].BigInt(&big.Int{}))
		gx.AddAssign(hy)
		parties[i] = Party{
			sKey1:   sKeys1[i],
			sKey2:   sKeys2[i],
			pKeyAff: *new(bls.G1Affine).FromJacobian(gx),
		}
	}

	// 计算价值向量承诺
	weightsF := make([]fr.Element, w.d)
	for i := 0; i < w.d; i++ {
		weightsF[i] = fr.NewElement(uint64(w.weights[i]))
	}
	w1Tau, _ := new(bls.G1Jac).MultiExp(w.PP.lagHTaus, weightsF, ecc.MultiExpConfig{})
	w2Tau, _ := new(bls.G2Jac).MultiExp(w.PP.lag2HTaus, weightsF, ecc.MultiExpConfig{})

	// 计算消除多项式承诺
	one := fr.One()
	var tauN fr.Element
	tauN.Exp(w.PP.tau, big.NewInt(int64(w.d)))
	tauN.Sub(&tauN, &one)
	vHTau := new(bls.G2Jac).ScalarMultiplication(&w.PP.g2, tauN.BigInt(&big.Int{}))

	// 计算g^{1/d},h^{1/d}
	nInv := fr.NewElement(uint64(w.d))
	nInv.Inverse(&nInv)
	gNInv := *new(bls.G2Affine).FromJacobian(new(bls.G2Jac).ScalarMultiplication(&w.PP.g2, nInv.BigInt(&big.Int{})))
	hNInv := *new(bls.G2Affine).ScalarMultiplication(&w.PP.h2a, nInv.BigInt(&big.Int{}))

	w.vk = VK{
		w1Tau:  *new(bls.G1Affine).FromJacobian(w1Tau),
		w2Tau:  *new(bls.G2Affine).FromJacobian(w2Tau),
		vHTau:  *new(bls.G2Affine).FromJacobian(vHTau),
		g2NInv: gNInv,
		h2NInv: hNInv,
	}
	w.signers = parties
}

func (w *WTS) profitPf(Contribution []int) (bls.G1Affine, bls.G1Affine, bls.G1Affine) {
	//初始化三个向量
	bF := make([]fr.Element, w.d) // 贡献度向量
	wF := make([]fr.Element, w.d) // 价值向量
	rF := make([]fr.Element, w.d) // 内积向量
	for i := 0; i < w.d; i++ {
		bF[i] = fr.NewElement(uint64(Contribution[i]))
		wF[i] = fr.NewElement(uint64(w.weights[i]))
		rF[i].Mul(&bF[i], &wF[i])
	}

	//将原始点值域数据转换为多项式系数域
	w.PP.domain.FFTInverse(bF, fft.DIF)
	w.PP.domain.FFTInverse(wF, fft.DIF)
	w.PP.domain.FFTInverse(rF, fft.DIF)

	//计算系数域乘积
	w.PP.domain.FFT(bF, fft.DIT, true)
	w.PP.domain.FFT(wF, fft.DIT, true)
	w.PP.domain.FFT(rF, fft.DIT, true)

	//构造多项式差值部分
	one := fr.One()
	var den fr.Element
	den.Exp(w.PP.domain.FrMultiplicativeGen, big.NewInt(int64(w.PP.domain.Cardinality)))
	den.Sub(&den, &one).Inverse(&den)

	//构造差值多项式
	for i := 0; i < w.d; i++ {
		bF[i].Mul(&bF[i], &wF[i]).
			Sub(&bF[i], &rF[i]).
			Mul(&bF[i], &den)
	}

	//回到系数域
	w.PP.domain.FFTInverse(bF, fft.DIF, true)
	w.PP.domain.FFTInverse(rF, fft.DIF, true)
	fft.BitReverse(bF)
	fft.BitReverse(rF)

	//多项式指数映射
	qTau, _ := new(bls.G1Jac).MultiExp(w.PP.PoT, bF, ecc.MultiExpConfig{})
	rTau, _ := new(bls.G1Jac).MultiExp(w.PP.PoT[:w.d-1], rF[1:], ecc.MultiExpConfig{})
	pTauH, _ := new(bls.G1Jac).MultiExp(w.PP.PoTH, rF, ecc.MultiExpConfig{})

	qTauAff := new(bls.G1Affine).FromJacobian(qTau)
	rTauAff := new(bls.G1Affine).FromJacobian(rTau)
	pTauHAff := new(bls.G1Affine).FromJacobian(pTauH)

	return *qTauAff, *rTauAff, *pTauHAff
}

func (w *WTS) Sign(fid Message, m [][]fr.Element, typeVec []Message) (Tag, []bls.G1Affine, []bls.G1Affine) {
	n := len(m)

	// 初始化结果容器
	tag := Tag{
		t:     make([]bls.G1Affine, n),
		sigma: make([]bls.G1Affine, n),
		K:     make([]bls.G1Affine, n),
		xba:   make([]fr.Element, n),
		yba:   make([]fr.Element, n),
	}
	cTau := make([]bls.G1Affine, w.s)
	tdba := make([]bls.G1Affine, w.s)

	// 初始化用户贡献度向量
	c := make([][]int, w.s)
	for i := 0; i < w.s; i++ {
		c[i] = make([]int, w.d)
	}

	//为每一个数据m生成tag
	eta, _ := bls.HashToG1(fid, []byte{})
	var etaJac bls.G1Jac
	etaJac.FromAffine(&eta)
	for i := 0; i < n; i++ {
		//随机选择本次签名的用户和文件类型
		π := rand.Intn(w.s)
		mπ := m[i][π]
		l := rand.Intn(len(typeVec))
		c[π][l]++
		typeByte := typeVec[l]
		var rx, ry, hrx fr.Element
		rx.SetRandom()
		ry.SetRandom()
		rxbyte := rx.Bytes()
		hrx = HashTofr(rxbyte[:])
		fidi := concatMessageAndIndex(fid, i)
		hfi, _ := bls.HashToG1(fidi, []byte{})
		hfiJac := *new(bls.G1Jac).FromAffine(&hfi)

		//计算t_i
		htype := HashTofr(typeByte)
		exp := *new(fr.Element).Mul((&w.signers[π].sKey1), &htype)
		var tJac bls.G1Jac
		tJac.ScalarMultiplication(&etaJac, exp.BigInt(new(big.Int)))
		tag.t[i] = *new(bls.G1Affine).FromJacobian(&tJac)

		//计算sigma_i
		//第一部分
		one := fr.One()
		tmp1 := *new(fr.Element).Sub(&one, &htype)
		tmp1.Mul(&tmp1, &mπ).Mul(&tmp1, &w.signers[π].sKey1)
		expG := *new(fr.Element).Add(&rx, &tmp1)
		gTerm := *new(bls.G1Jac).ScalarMultiplication(&w.PP.g1, expG.BigInt(new(big.Int)))

		//第二部分
		hTerm := *new(bls.G1Jac).ScalarMultiplication(&w.PP.h1, ry.BigInt(new(big.Int)))

		//第三部分
		exph := *new(fr.Element).Mul(&htype, &hrx)
		pkπJac := *new(bls.G1Jac).FromAffine(&w.signers[π].pKeyAff)
		pkπTerm := *new(bls.G1Jac).ScalarMultiplication(&pkπJac, exph.BigInt(new(big.Int)))

		//第四部分
		var prodTerm bls.G1Jac
		for j := 0; j < w.s; j++ {
			if j == π {
				continue
			}
			tmp2 := *new(bls.G1Jac).FromAffine(&w.signers[j].pKeyAff)
			tmp2.ScalarMultiplication(&tmp2, m[i][j].BigInt(new(big.Int)))
			prodTerm.AddAssign(&tmp2)
		}

		// 合并所有部分
		var sigma bls.G1Jac
		sigma.Set(&gTerm)
		sigma.AddAssign(&hTerm)
		sigma.AddAssign(&hfiJac)
		sigma.AddAssign(&pkπTerm)
		sigma.AddAssign(&prodTerm)
		tag.sigma[i] = *new(bls.G1Affine).FromJacobian(&sigma)

		//计算K
		var er bls.G1Jac
		er.ScalarMultiplication(&etaJac, rx.BigInt(new(big.Int)))
		var sum fr.Element
		for j := 0; j < w.s; j++ {
			if j == π {
				continue
			}
			sum.Add(&sum, &m[i][j])
		}
		expT := *new(fr.Element).Add(&hrx, &sum)
		var thm bls.G1Jac
		thm.ScalarMultiplication(&tJac, expT.BigInt(new(big.Int)))
		er.AddAssign(&thm)
		er.AddAssign(&hfiJac)
		tag.K[i] = *new(bls.G1Affine).FromJacobian(&er)

		//计算xba,yba
		var mr fr.Element
		mr.Sub(&mπ, &hrx)
		var xhmh fr.Element
		xhmh.Mul(&w.signers[π].sKey1, &htype).Mul(&xhmh, &mr)
		var xba fr.Element
		xba.Sub(&rx, &xhmh)
		tag.xba[i] = xba

		var tr fr.Element
		tr.Mul(&htype, &hrx)
		var mtr fr.Element
		mtr.Sub(&mπ, &tr)
		var ymtr fr.Element
		ymtr.Mul(&w.signers[π].sKey2, &mtr)
		var yba fr.Element
		yba.Sub(&ry, &ymtr)
		tag.yba[i] = yba
	}

	for i := 0; i < w.s; i++ {
		//计算用户权益和贡献度承诺
		W := 0
		contributionF := make([]fr.Element, w.d)
		for j := 0; j < w.d; j++ {
			W += c[i][j] * w.weights[j]
			contributionF[j] = fr.NewElement(uint64(c[i][j]))
		}
		cTaui, _ := new(bls.G1Jac).MultiExp(w.PP.lagHTaus, contributionF, ecc.MultiExpConfig{})
		cTau[i] = *new(bls.G1Affine).FromJacobian(cTaui)

		//计算加密链接陷门
		var td bls.G1Jac
		td.ScalarMultiplication(&etaJac, w.signers[i].sKey1.BigInt(new(big.Int)))
		tdba[i] = *new(bls.G1Affine).FromJacobian(&td)
	}
	return tag, cTau, tdba
}

func (w *WTS) ChalGen(c int) Chal {
	// 初始化挑战容器
	chal := Chal{
		index: make([]int, c),
		v:     make([]fr.Element, c),
	}
	chal.index = getRandomIndices(w.n, c)

	//v
	for i := 0; i < c; i++ {
		chal.v[i].SetRandom()
	}
	return chal
}

func (w *WTS) ProofGen(c int, m [][]fr.Element, tag Tag, chal Chal) Proof {
	//计算证据
	//初始化证据容器
	proof := Proof{
		Psigma: bls.G1Affine{},
		X:      fr.Element{},
		Y:      fr.Element{},
		mu:     make([]fr.Element, w.s),
		theta:  bls.G1Affine{},
	}

	PsigmaProd := *new(bls.G1Jac)
	var X, Y fr.Element
	muProd := make([]fr.Element, w.s)
	theta := *new(bls.G1Jac)

	for i := 0; i < c; i++ {
		// 获取挑战索引
		index := chal.index[i]
		vi := chal.v[i]

		//Psigma
		sigma_i := tag.sigma[index]
		K_i := tag.K[index]
		sigmaJac := *new(bls.G1Jac).FromAffine(&sigma_i)
		kJac := *new(bls.G1Jac).FromAffine(&K_i)

		PsigmaTerm := sigmaJac
		PsigmaTerm.AddAssign(&kJac)
		PsigmaTerm.ScalarMultiplication(&PsigmaTerm, vi.BigInt(new(big.Int)))
		PsigmaProd.AddAssign(&PsigmaTerm)

		//X
		var xTerm fr.Element
		xTerm.Mul(&vi, &tag.xba[index])
		X.Add(&X, &xTerm)

		//Y
		var yTerm fr.Element
		yTerm.Mul(&vi, &tag.yba[index])
		Y.Add(&Y, &yTerm)

		//mu,theta
		var thetaProd fr.Element
		var thetaJac bls.G1Jac
		muTerm := make([]fr.Element, w.s)
		for j := 0; j < w.s; j++ {
			muTerm[j].Mul(&vi, &m[index][j])
			muProd[j].Add(&muProd[j], &muTerm[j])
			thetaProd.Add(&thetaProd, &muTerm[j])
		}
		ti := *new(bls.G1Jac).FromAffine(&tag.t[index])
		thetaJac.ScalarMultiplication(&ti, thetaProd.BigInt(new(big.Int)))
		theta.AddAssign(&thetaJac)
	}

	proof.Psigma = *new(bls.G1Affine).FromJacobian(&PsigmaProd)
	proof.X = X
	proof.Y = Y
	proof.mu = muProd
	proof.theta = *new(bls.G1Affine).FromJacobian(&theta)

	return proof
}

func (w *WTS) Verify(chal Chal, proof Proof, fid Message) bool {

	//等式1左边部分
	var left bls.G1Jac
	eta, _ := bls.HashToG1(fid, []byte{})
	etaJac := *new(bls.G1Jac).FromAffine(&eta)
	getaX := w.PP.g1
	getaX.AddAssign(&etaJac)
	getaX.ScalarMultiplication(&getaX, proof.X.BigInt(new(big.Int)))
	left.AddAssign(&getaX)

	var hY bls.G1Jac
	hY.ScalarMultiplication(&w.PP.h1, proof.Y.BigInt(new(big.Int)))
	left.AddAssign(&hY)

	Theta := *new(bls.G1Jac).FromAffine(&proof.theta)
	left.AddAssign(&Theta)

	c := len(chal.index)
	var hfiterm bls.G1Jac
	for i := 0; i < c; i++ {
		index := chal.index[i]
		vi := chal.v[i]
		fidi := concatMessageAndIndex(fid, index)
		hfi, _ := bls.HashToG1(fidi, []byte{})
		hfiJac := *new(bls.G1Jac).FromAffine(&hfi)
		hfiJac.ScalarMultiplication(&hfiJac, vi.BigInt(new(big.Int)))
		hfiterm.AddAssign(&hfiJac)
	}
	left.AddAssign(&hfiterm)
	left.AddAssign(&hfiterm)

	for i := 0; i < w.s; i++ {
		var pkmu bls.G1Jac
		pkJac := *new(bls.G1Jac)
		pkJac.FromAffine(&w.signers[i].pKeyAff)
		pkmu.ScalarMultiplication(&pkJac, proof.mu[i].BigInt(new(big.Int)))
		left.AddAssign(&pkmu)
	}
	leftAFF := *new(bls.G1Affine).FromJacobian(&left)

	// 验证等式
	if leftAFF.Equal(&proof.Psigma) {
		return true
	} else {
		return false
	}

}

func (w *WTS) ProfitGen(ti []bls.G1Affine, td bls.G1Affine, typeVec []Message) Profit {
	//初始化结构体
	profit := Profit{
		cTau: bls.G1Affine{},
		P:    0,
		pi: IPAProof{
			qTau: bls.G1Affine{},
			rTau: bls.G1Affine{},
			pTau: bls.G1Affine{},
		},
	}

	c := make([]int, w.d)

	//构造映射表
	targetMap := make(map[[48]byte]int)
	for i := 0; i < w.d; i++ {
		tdJac := *new(bls.G1Jac).FromAffine(&td)
		htype := HashTofr(typeVec[i])
		tdJac.ScalarMultiplication(&tdJac, htype.BigInt(new(big.Int)))
		tdAff := *new(bls.G1Affine).FromJacobian(&tdJac)
		key := tdAff.Bytes()
		targetMap[key] = i
	}

	//链接生成用户贡献向量
	for i := 0; i < w.n; i++ {
		key := ti[i].Bytes()
		if j, ok := targetMap[key]; ok {
			c[j]++
		}
	}

	//计算用户权益
	W := 0
	for i := 0; i < w.d; i++ {
		W += c[i] * w.weights[i]
	}
	profit.P = W

	//计算内积论证证据
	var qwTau, rwTau, pwTauH bls.G1Affine
	qwTau, rwTau, pwTauH = w.profitPf(c)
	profit.pi.qTau = qwTau
	profit.pi.rTau = rwTau
	profit.pi.pTau = pwTauH

	return profit
}

func (w *WTS) ProfitVerify(profit Profit) bool {
	//验证内积结果的正确性
	W := fr.NewElement(uint64(profit.P))
	gw := new(bls.G1Affine).ScalarMultiplication(&w.PP.g1a, W.BigInt(&big.Int{}))

	// 验证配对1
	lhs1, _ := bls.Pair([]bls.G1Affine{profit.cTau}, []bls.G2Affine{w.vk.w2Tau})
	rhs1, _ := bls.Pair([]bls.G1Affine{profit.pi.qTau, profit.pi.rTau, *gw}, []bls.G2Affine{w.vk.vHTau, w.PP.g2Tau, w.vk.g2NInv})
	res1 := lhs1.Equal(&rhs1)
	if res1 {
		/*fmt.Println("配对1验证通过")*/
	} else {
		/*fmt.Println("配对1验证失败")*/
		return false
	}

	// 验证配对2
	lhs2, _ := bls.Pair([]bls.G1Affine{profit.pi.pTau}, []bls.G2Affine{w.PP.g2a})
	rhs2, _ := bls.Pair([]bls.G1Affine{profit.pi.rTau, *gw}, []bls.G2Affine{w.PP.hTauHAff, w.vk.h2NInv})
	res2 := lhs2.Equal(&rhs2)
	if res2 {
		/*fmt.Println("配对2验证通过")*/
		return true
	} else {
		/*fmt.Println("配对2验证失败")*/
		return false
	}
}

// 获得系数
func (w *WTS) getFSChal(vals []bls.G1Affine, ths int) fr.Element {
	n := len(vals)
	hMsg := make([]byte, n*48+4)
	for i, val := range vals {
		mBytes := val.Bytes()
		copy(hMsg[i*48:(i+1)*48], mBytes[:])
	}
	tBytes := make([]byte, 4)
	binary.LittleEndian.PutUint32(tBytes, uint32(ths))
	copy(hMsg[n*48:], tBytes)

	hFunc := sha256.New()
	hFunc.Reset()
	return *new(fr.Element).SetBytes(hFunc.Sum(hMsg))
}

func (w *WTS) BatchProfitGen(ti []bls.G1Affine, typeVec []Message, td []bls.G1Affine) BatchProfit {
	//初始化结构体
	profit := BatchProfit{
		cTau: make([]bls.G1Affine, w.s),
		P:    make([]int, w.s),
		pi: IPAProof{
			qTau: bls.G1Affine{},
			rTau: bls.G1Affine{},
			pTau: bls.G1Affine{},
		},
	}

	// 计算所有用户的贡献向量
	c := make([][]int, w.s)
	for i := 0; i < w.s; i++ {
		c[i] = make([]int, w.d)
	}

	//计算用户链接映射表
	type key struct{ i, j int }
	targetMap := make(map[[48]byte]key)
	for i := 0; i < w.s; i++ {
		for j := 0; j < w.d; j++ {
			tdJac := *new(bls.G1Jac).FromAffine(&td[i])
			htype := HashTofr(typeVec[j])
			tdJac.ScalarMultiplication(&tdJac, htype.BigInt(new(big.Int)))
			tdAff := *new(bls.G1Affine).FromJacobian(&tdJac)
			targetMap[tdAff.Bytes()] = key{i, j}
		}
	}

	//链接计算所有用户的贡献向量
	for i := 0; i < w.n; i++ {
		k := ti[i].Bytes()
		if ij, ok := targetMap[k]; ok {
			c[ij.i][ij.j]++
		}
	}

	//计算并聚合证据
	var qoTau, roTau, poTauH bls.G1Jac
	for i := 0; i < w.s; i++ {
		//计算用户权益
		W := 0
		contributionF := make([]fr.Element, w.d)
		for j := 0; j < w.d; j++ {
			W += c[i][j] * w.weights[j]
			contributionF[j] = fr.NewElement(uint64(c[i][j]))
		}
		cTaui, _ := new(bls.G1Jac).MultiExp(w.PP.lagHTaus, contributionF, ecc.MultiExpConfig{})
		cTau := *new(bls.G1Affine).FromJacobian(cTaui)
		profit.P[i] = W

		//聚合证据
		var qwTau, rwTau, pwTauH bls.G1Affine
		qwTau, rwTau, pwTauH = w.profitPf(c[i])
		qwTauJac := *new(bls.G1Jac).FromAffine(&qwTau)
		rwTauJac := *new(bls.G1Jac).FromAffine(&rwTau)
		pwTauHJac := *new(bls.G1Jac).FromAffine(&pwTauH)

		xi := w.getFSChal([]bls.G1Affine{cTau, w.vk.w1Tau}, W)
		xiInt := xi.BigInt(&big.Int{})

		qoTau.AddAssign(qwTauJac.ScalarMultiplication(&qwTauJac, xiInt))
		roTau.AddAssign(rwTauJac.ScalarMultiplication(&rwTauJac, xiInt))
		poTauH.AddAssign(pwTauHJac.ScalarMultiplication(&pwTauHJac, xiInt))
	}

	profit.pi.qTau = *new(bls.G1Affine).FromJacobian(&qoTau)
	profit.pi.rTau = *new(bls.G1Affine).FromJacobian(&roTau)
	profit.pi.pTau = *new(bls.G1Affine).FromJacobian(&poTauH)

	return profit
}

func (w *WTS) BatchProfitVerify(profit BatchProfit) bool {

	//  获取系数并聚合
	var oTau, gw bls.G1Jac
	for i := 0; i < w.s; i++ {
		var oTauTerm, gwTerm bls.G1Jac
		xi := w.getFSChal([]bls.G1Affine{profit.cTau[i], w.vk.w1Tau}, profit.P[i])
		cTau := *new(bls.G1Jac).FromAffine(&profit.cTau[i])
		oTauTerm.ScalarMultiplication(&cTau, xi.BigInt(&big.Int{}))
		oTau.AddAssign(&oTauTerm)

		Pfr := fr.NewElement(uint64(profit.P[i]))
		xiP := *new(fr.Element).Mul(&xi, &Pfr)
		gwTerm.ScalarMultiplication(&w.PP.g1, xiP.BigInt(&big.Int{}))
		gw.AddAssign(&gwTerm)
	}

	oTauAff := *new(bls.G1Affine).FromJacobian(&oTau)
	gwAFF := new(bls.G1Affine).FromJacobian(&gw)

	//批量验证配对1
	lhs1, _ := bls.Pair([]bls.G1Affine{oTauAff}, []bls.G2Affine{w.vk.w2Tau})
	rhs1, _ := bls.Pair([]bls.G1Affine{profit.pi.qTau, profit.pi.rTau, *gwAFF}, []bls.G2Affine{w.vk.vHTau, w.PP.g2Tau, w.vk.g2NInv})
	res1 := lhs1.Equal(&rhs1)
	if res1 {
		/*fmt.Println("批量配对1验证通过")*/
	} else {
		/*fmt.Println("批量配对1验证失败")*/
		return false
	}

	//批量验证配对2
	lhs2, _ := bls.Pair([]bls.G1Affine{profit.pi.pTau}, []bls.G2Affine{w.PP.g2a})
	rhs2, _ := bls.Pair([]bls.G1Affine{profit.pi.rTau, *gwAFF}, []bls.G2Affine{w.PP.hTauHAff, w.vk.h2NInv})
	res2 := lhs2.Equal(&rhs2)
	if res2 {
		/*fmt.Println("批量配对2验证通过")*/
		return true
	} else {
		/*fmt.Println("批量配对2验证失败")*/
		return false
	}
}

func (w *WTS) QueryGen(Contribution []int, ell int) (int, bls.G1Affine) {
	psi := make([]fr.Element, w.d)
	cjl := Contribution[ell]
	cjlfr := fr.NewElement(uint64(cjl))
	ellfr := fr.NewElement(uint64(ell))
	for i := 0; i < w.d; i++ {
		psi[i] = fr.NewElement(uint64(Contribution[i]))
	}
	w.PP.domain.FFTInverse(psi, fft.DIF)
	w.PP.domain.FFT(psi, fft.DIT, true)
	for i := 0; i < w.d; i++ {
		psi[i].Sub(&psi[i], &cjlfr).Div(&psi[i], &ellfr)
	}
	w.PP.domain.FFTInverse(psi, fft.DIF, true)
	fft.BitReverse(psi)
	piell, _ := new(bls.G1Jac).MultiExp(w.PP.PoT, psi, ecc.MultiExpConfig{})
	piellAff := new(bls.G1Affine).FromJacobian(piell)
	return cjl, *piellAff
}

func (w *WTS) QueryVerify(gcj bls.G1Affine, pi bls.G1Affine, ell int, cjl int) bool {
	lhs, _ := bls.Pair([]bls.G1Affine{gcj}, []bls.G2Affine{w.PP.g2a})

	ellfr := fr.NewElement(uint64(ell))
	gell := *new(bls.G2Jac).ScalarMultiplication(&w.PP.g2, ellfr.BigInt(&big.Int{}))
	gtauell := *new(bls.G2Jac).FromAffine(&w.PP.g2Tau)
	gtauell.SubAssign(&gell)
	gtlAff := *new(bls.G2Affine).FromJacobian(&gtauell)
	rhs_1, _ := bls.Pair([]bls.G1Affine{pi}, []bls.G2Affine{gtlAff})

	rhs_2, _ := bls.Pair([]bls.G1Affine{w.PP.g1a}, []bls.G2Affine{w.PP.g2a})
	rhs_2.Exp(rhs_1, ellfr.BigInt(&big.Int{}))
	rhs_1.Mul(&rhs_1, &rhs_2)

	res := lhs.Equal(&rhs_1)
	if res {
		return true
	} else {
		return false
	}
}
