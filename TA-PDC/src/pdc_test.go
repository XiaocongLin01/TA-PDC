package pdc

import (
	"flag"
	"fmt"
	"math/big"
	"testing"

	"github.com/consensys/gnark-crypto/ecc"
	bls "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
)

var NUM_NODES = flag.Int("type", 1<<4, "Number of Types")

// 用户密钥生成算法正确性测试
func TestKeyGen(t *testing.T) {
	s := 4
	d := 4
	n := 128

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)
	w := NewWTS(d, s, n, weights, PP)

	//查看并检查用户公私钥对
	for i, party := range w.signers {
		xBig := party.sKey1.BigInt(&big.Int{})
		yBig := party.sKey2.BigInt(&big.Int{})

		var g1Jac, h1Jac bls.G1Jac
		g1Jac.ScalarMultiplication(&w.PP.g1, xBig)
		h1Jac.ScalarMultiplication(&w.PP.h1, yBig)

		var pkRecomputed bls.G1Jac
		pkRecomputed.Set(&g1Jac)
		pkRecomputed.AddAssign(&h1Jac)

		var pkAffine bls.G1Affine
		pkAffine.FromJacobian(&pkRecomputed)
		if !pkAffine.Equal(&party.pKeyAff) {
			t.Errorf("Public key mismatch at index %d", i)
		}
	}

}

// 审计算法正确性测试
func TestAudit(t *testing.T) {
	s := 8
	d := 8
	n := 1024

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)
	w := NewWTS(d, s, n, weights, PP)

	//构造 m 二维数组
	m := make([][]fr.Element, n)
	for i := 0; i < n; i++ {
		m[i] = make([]fr.Element, s)
		for j := 0; j < s; j++ {
			m[i][j].SetUint64(uint64((i + j + 1) % 5))
		}
	}

	typeVec := []Message{
		[]byte("image"),
		[]byte("video"),
		[]byte("txt"),
		[]byte("pdf"),
	}

	// 构造 fid 和 event
	fid := []byte("file-1")

	tag, _, _ := w.Sign(fid, m, typeVec)

	chal := w.ChalGen(32)

	proof := w.ProofGen(32, m, tag, chal)

	result := w.Verify(chal, proof, fid)
	if result {
		fmt.Println("验证通过")
	} else {
		fmt.Println("验证失败")
	}
}

// 单个用户权益正确性测试
func TestProfit(t *testing.T) {
	s := 16
	d := 8
	n := 4096

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)
	w := NewWTS(d, s, n, weights, PP)

	//构造 m 二维数组
	m := make([][]fr.Element, n)
	for i := 0; i < n; i++ {
		m[i] = make([]fr.Element, s)
		for j := 0; j < s; j++ {
			m[i][j].SetUint64(uint64((i + j + 1) % 5))
		}
	}

	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}

	// 构造 fid
	fid := []byte("file-1")

	//贡献向量
	var ctau []bls.G1Affine
	var tag Tag
	var td []bls.G1Affine
	tag, ctau, td = w.Sign(fid, m, typeVec)

	profit := w.ProfitGen(tag.t, td[2], typeVec)
	profit.cTau = ctau[2]
	result := w.ProfitVerify(profit)
	if result {
		fmt.Println("配对通过")
	} else {
		fmt.Println("配对失败")
	}
}

// 批量权益正确性测试
func TestBatchProfit(t *testing.T) {
	s := 16
	d := 8
	n := 4096

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)
	w := NewWTS(d, s, n, weights, PP)

	//构造 m 二维数组
	m := make([][]fr.Element, n)
	for i := 0; i < n; i++ {
		m[i] = make([]fr.Element, s)
		for j := 0; j < s; j++ {
			m[i][j].SetUint64(uint64((i + j + 1) % 5))
		}
	}

	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}

	// 构造 fid
	fid := []byte("file-1")

	//贡献向量
	var ctau []bls.G1Affine
	var tag Tag
	var td []bls.G1Affine
	tag, ctau, td = w.Sign(fid, m, typeVec)

	profit := w.BatchProfitGen(tag.t, typeVec, td)
	profit.cTau = ctau
	result := w.BatchProfitVerify(profit)

	if result {
		fmt.Println("批量配对通过")
	} else {
		fmt.Println("批量配对失败")
	}
}

// 获得群元素大小和各种运算时间
func BenchmarkGroupOperations(b *testing.B) {
	// 初始化 fr 元素
	var frElement fr.Element
	frElement.SetRandom()

	// 初始化 G1 和 G2 元素
	g1Jac, g2Jac, g1Aff, g2Aff := bls.Generators()

	fmt.Printf("fr 元素比特长度: %d bits\n", frElement.BitLen())
	fmt.Printf("G1 元素比特长度: %d bits\n", g1Aff.X.BitLen())

	//测试G1乘法时间
	b.Run("G1Mul Time", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			g1Jac.AddAssign(&g1Jac)
		}
	})

	//测试G2乘法时间
	b.Run("G2Mul Time", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			g2Jac.AddAssign(&g2Jac)
		}
	})

	//测试G1幂乘时间
	b.Run("G1Exp Time", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			g1Jac.ScalarMultiplication(&g1Jac, frElement.BigInt(&big.Int{}))
		}
	})

	//测试G2幂乘时间
	b.Run("G2Exp Time", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			g2Jac.ScalarMultiplication(&g2Jac, frElement.BigInt(&big.Int{}))
		}
	})

	//测试配对时间
	b.Run("Pairing Time", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			bls.Pair([]bls.G1Affine{g1Aff}, []bls.G2Affine{g2Aff})
		}
	})

	//测试GT乘法时间
	b.Run("GTMul Time", func(b *testing.B) {
		gt, _ := bls.Pair([]bls.G1Affine{g1Aff}, []bls.G2Affine{g2Aff})
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			gt.Mul(&gt, &gt)
		}
	})

	//测试G1Hash时间
	b.Run("G1Hash Time", func(b *testing.B) {
		str := []byte("test")
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			bls.HashToG1(str, []byte{})
		}
	})
}

// 测试Setup阶段时间
func BenchmarkSetup(b *testing.B) {
	d := []int{4, 8, 12, 16, 20}

	for _, value := range d {
		b.Run("Setup Time - "+fmt.Sprintf("%d", value)+" value type", func(b *testing.B) {
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				Setup(value)
			}
		})
	}
}

// 测试KeyGen阶段时间
func BenchmarkKenGen(b *testing.B) {
	d := 8
	s := []int{4, 8, 12, 16, 20}
	n := 1024

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)

	for _, value := range s {
		b.Run("KeyGen Time - "+fmt.Sprintf("%d", value)+" users", func(b *testing.B) {
			w := WTS{
				d:       d,
				n:       n,
				s:       value,
				weights: weights,
				PP:      PP,
			}
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.keyGen()
			}
		})
	}
}

// 测试Sign阶段时间
func BenchmarkSign(b *testing.B) {
	d := 8
	/*s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}*/
	s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}
	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}
	PP := Setup(d)
	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}
	fid := []byte("file-1")

	for _, value := range n {
		w := NewWTS(d, s, value, weights, PP)
		m := make([][]fr.Element, value)
		for i := 0; i < value; i++ {
			m[i] = make([]fr.Element, s)
			for j := 0; j < s; j++ {
				m[i][j].SetUint64(uint64((i + j + 1) % 5))
			}
		}
		b.Run("Sign Time - "+fmt.Sprintf("%d", ((value*s*32)/1024))+" KB", func(b *testing.B) {
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.Sign(fid, m, typeVec)
			}
		})
	}

}

// 测试证据生成时间
func BenchmarkProofGen(b *testing.B) {
	d := 8
	/*s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}*/
	s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}
	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}
	PP := Setup(d)
	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}
	fid := []byte("file-1")

	for _, value := range n {
		b.Run("ProofGen Time - "+fmt.Sprintf("%d", ((value*s*32)/1024))+" KB", func(b *testing.B) {
			w := NewWTS(d, s, value, weights, PP)
			m := make([][]fr.Element, value)
			for i := 0; i < value; i++ {
				m[i] = make([]fr.Element, s)
				for j := 0; j < s; j++ {
					m[i][j].SetUint64(uint64((i + j + 1) % 5))
				}
			}

			tag, _, _ := w.Sign(fid, m, typeVec)
			chal := w.ChalGen(460)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.ProofGen(460, m, tag, chal)
			}
		})
	}

}

// 测试审计验证时间
func BenchmarkProofVerify(b *testing.B) {
	d := 8
	s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}
	/*s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}*/
	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}
	PP := Setup(d)
	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}
	fid := []byte("file-1")

	for _, value := range n {
		b.Run("ProofVerify Time - "+fmt.Sprintf("%d", ((value*s*32)/1024))+" KB", func(b *testing.B) {
			w := NewWTS(d, s, value, weights, PP)
			m := make([][]fr.Element, value)
			for i := 0; i < value; i++ {
				m[i] = make([]fr.Element, s)
				for j := 0; j < s; j++ {
					m[i][j].SetUint64(uint64((i + j + 1) % 5))
				}
			}

			tag, _, _ := w.Sign(fid, m, typeVec)
			chal := w.ChalGen(460)
			proof := w.ProofGen(460, m, tag, chal)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.Verify(chal, proof, fid)
			}
		})
	}

}

// 测试权益证据生成时间
func BenchmarkProfitGen(b *testing.B) {
	d := 500
	// s := []int{4, 8, 12, 16, 20, 24, 28, 32}
	s := []int{1000}
	n := 1024
	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)

	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}

	// 构造 fid 和 event
	fid := []byte("file-1")

	for _, value := range s {
		w := NewWTS(d, value, n, weights, PP)
		m := make([][]fr.Element, n)
		for i := 0; i < n; i++ {
			m[i] = make([]fr.Element, value)
			for j := 0; j < value; j++ {
				m[i][j].SetUint64(uint64((i + j + 1) % 5))
			}
		}

		var tag Tag
		var td []bls.G1Affine
		tag, _, td = w.Sign(fid, m, typeVec)

		b.Run("ProfitGen Time - "+fmt.Sprintf("%d", value)+" users", func(b *testing.B) {
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.BatchProfitGen(tag.t, typeVec, td)
			}
		})
	}

}

// 测试权益验证时间
func BenchmarkProfitVerify(b *testing.B) {
	d := 8
	s := []int{4, 8, 12, 16, 20, 24, 28, 32}
	// s := []int{1000}
	n := 1024
	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)

	var typeVec []Message
	for i := 0; i < d; i++ {
		typeVec = append(typeVec, Message([]byte("type"+fmt.Sprintf("%d", i))))
	}

	// 构造 fid 和 event
	fid := []byte("file-1")

	for _, value := range s {
		w := NewWTS(d, value, n, weights, PP)
		m := make([][]fr.Element, n)
		for i := 0; i < n; i++ {
			m[i] = make([]fr.Element, value)
			for j := 0; j < value; j++ {
				m[i][j].SetUint64(uint64((i + j + 1) % 5))
			}
		}

		var ctau []bls.G1Affine
		var tag Tag
		var td []bls.G1Affine
		tag, ctau, td = w.Sign(fid, m, typeVec)
		profit := w.BatchProfitGen(tag.t, typeVec, td)
		profit.cTau = ctau

		b.Run("ProfitVergy Time - "+fmt.Sprintf("%d", value)+" users", func(b *testing.B) {
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				w.BatchProfitVerify(profit)
			}
		})
	}

}

func BenchmarkQuery(b *testing.B) {
	d := 8
	s := 16
	n := 1024
	ell := 6
	contribution := []int{1, 2, 3, 4, 5, 6, 7, 8}

	weights := make([]int, d)
	for i := 0; i < d; i++ {
		weights[i] = i + 1
	}

	PP := Setup(d)

	w := NewWTS(d, s, n, weights, PP)

	psi := make([]fr.Element, w.d)
	for i := 0; i < w.d; i++ {
		psi[i] = fr.NewElement(uint64(contribution[i]))
	}
	cTau, _ := new(bls.G1Jac).MultiExp(w.PP.lagHTaus, psi, ecc.MultiExpConfig{})
	cTauAff := *new(bls.G1Affine).FromJacobian(cTau)

	b.Run("Query Time - "+fmt.Sprintf("%d", s)+" number", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			cjl, pi := w.QueryGen(contribution, ell)
			w.QueryVerify(cTauAff, pi, ell, cjl)
		}
	})

}
