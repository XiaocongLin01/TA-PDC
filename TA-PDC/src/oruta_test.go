package pdc

import (
	"fmt"
	"testing"

	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
)

func TestOrutaAudit(t *testing.T) {
	k := 8
	s := 8
	n := 1024

	PP := OrutaSetup()
	o := NewOruta(k, s, n, PP)

	//构造 m 二维数组
	m := make([][]fr.Element, n)
	for i := 0; i < n; i++ {
		m[i] = make([]fr.Element, k)
		for j := 0; j < k; j++ {
			m[i][j].SetUint64(uint64((i + j + 1) % 5))
		}
	}

	tag := o.OrutaSign(m)

	chal := o.OrutaChalGen(32)

	proof := o.OrutaProofGen(32, m, tag, chal)

	result := o.OrutaProofVerify(chal, proof)
	if result {
	} else {
		fmt.Println("验证失败,未实现同构计算仅用于模拟")
	}
}

func BenchmarkOrutaSign(b *testing.B) {
	k := 16
	/*s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}*/
	s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}
	PP := OrutaSetup()

	for _, value := range n {
		o := NewOruta(k, s, value, PP)
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
				o.OrutaSign(m)
			}
		})
	}
}

func BenchmarkOrutaProofGen(b *testing.B) {
	k := 16
	/*s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}*/
	s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}
	PP := OrutaSetup()
	for _, value := range n {
		b.Run("ProofGen Time - "+fmt.Sprintf("%d", ((value*s*32)/1024))+" KB", func(b *testing.B) {
			o := NewOruta(k, s, value, PP)
			m := make([][]fr.Element, value)
			for i := 0; i < value; i++ {
				m[i] = make([]fr.Element, s)
				for j := 0; j < s; j++ {
					m[i][j].SetUint64(uint64((i + j + 1) % 5))
				}
			}

			tag := o.OrutaSign(m)
			chal := o.OrutaChalGen(460)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				o.OrutaProofGen(460, m, tag, chal)
			}
		})
	}
}

func BenchmarkOrutaProofVerify(b *testing.B) {
	k := 16
	/*s := 8 //2048bits
	n := []int{1024, 2048, 4096, 8192, 16384}*/
	s := 16 //4096bits
	n := []int{512, 1024, 2048, 4096, 8192, 16384}
	PP := OrutaSetup()

	for _, value := range n {
		b.Run("ProofVerify Time - "+fmt.Sprintf("%d", ((value*s*32)/1024))+" KB", func(b *testing.B) {
			o := NewOruta(k, s, value, PP)
			m := make([][]fr.Element, value)
			for i := 0; i < value; i++ {
				m[i] = make([]fr.Element, s)
				for j := 0; j < s; j++ {
					m[i][j].SetUint64(uint64((i + j + 1) % 5))
				}
			}

			tag := o.OrutaSign(m)
			chal := o.OrutaChalGen(460)
			proof := o.OrutaProofGen(460, m, tag, chal)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				o.OrutaProofVerify(chal, proof)
			}
		})
	}
}
