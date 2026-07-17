# context.md — 2026 암호분석경진대회 7번 (RSA 부분키 노출 → 격자 공격)

이 문서는 **주어진 p의 일부 비트 + 공개키(N,e) + 암호문 ct** 로부터 N을 소인수분해하고
교과서적 RSA를 복호화하는 코드를 바로 작성·실행하기 위한 레퍼런스다.

---

## 1. 문제 요약

- N = pq, p·q 각각 1024비트, N은 2048비트.
- p의 일부 비트가 유출됨. `MASK`의 비트가 1인 위치가 "아는 비트"이고, 그 값은 `p & MASK`.
- 공격: 미지 비트를 격자(Herrmann–May / 다변수 Coppersmith)로 복원 → p 획득 → 표준 RSA 복호화.

핵심 수치:
- 아는 비트 672개, 모르는 비트 352개.
- 모르는 비트는 7개 블록. 그중 4비트짜리 2개는 **전수조사(brute force)**, 나머지 5개 블록(344비트)은 격자로 복원.
- 약수 크기 지수 **β = 1024/2048 = 0.5**. n=5 블록에서 Herrmann–May 점근 허용치 ≈ 361비트 > 344비트 → 복원 가능.

---

## 2. 주어진 데이터 (그대로 복붙해서 사용)

```python
N = int("""
e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1
""".replace("\n","").replace(" ",""), 16)

e = 0x10001

ct = int("""
8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b
""".replace("\n","").replace(" ",""), 16)

MASK = int("""
ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff
""".replace("\n","").replace(" ",""), 16)

P_AND_MASK = int("""
ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9
""".replace("\n","").replace(" ",""), 16)
```

> 주의: 위 4개 값은 모두 1024/2048비트 정수. 비트 인덱스는 LSB=0 관례를 쓴다.
> (`(x >> i) & 1` 이 bit i). MASK/p&MASK는 1024비트, N/ct는 2048비트.

---

## 3. 미지 블록 (MASK 파싱 결과)

MASK에서 0인 연속 구간 = 모르는 비트 블록. 파싱 코드:

```python
bits = [(MASK >> i) & 1 for i in range(1024)]
runs = []
i = 0
while i < 1024:
    if bits[i] == 0:
        j = i
        while j < 1024 and bits[j] == 0:
            j += 1
        runs.append((i, j - i))   # (offset, length)
        i = j
    else:
        i += 1
# runs = [(offset, length), ...]
```

결과 (offset = 블록의 최하위 비트 위치, length = 비트 수):

| # | offset | length | 처리 |
|---|--------|--------|------|
| 0 | 150 | 4  | **brute force** |
| 1 | 265 | 84 | 격자 변수 x1 |
| 2 | 362 | 58 | 격자 변수 x2 |
| 3 | 600 | 69 | 격자 변수 x3 |
| 4 | 682 | 87 | 격자 변수 x4 |
| 5 | 784 | 46 | 격자 변수 x5 |
| 6 | 920 | 4  | **brute force** |

- 격자 변수 5개, 총 344비트. 각 변수 bound `X_i = 2^length`.
- 4비트 블록 2개(offset 150, 920)는 256가지 전수조사.

---

## 4. 공격 모델 (Herrmann–May)

p를 미지 블록의 선형식으로 표현:

```
p = p_base + b150·2^150 + x1·2^265 + x2·2^362 + x3·2^600 + x4·2^682 + x5·2^784 + b920·2^920
```

- `p_base` = P_AND_MASK (아는 비트만 채워진 값; 미지 자리는 0). 실제로 P_AND_MASK가 곧 p_base다(미지 자리는 이미 0).
- `b150 ∈ [0,16)`, `b920 ∈ [0,16)` = 전수조사 대상.
- 전수조사로 b150,b920을 고정하면 상수항
  `c = P_AND_MASK + b150·2^150 + b920·2^920` 가 정해지고,
  남은 식은 5변수 선형:

```
f(x1..x5) = c + x1·2^265 + x2·2^362 + x3·2^600 + x4·2^682 + x5·2^784
```

이 f는 참값에서 `f ≡ 0 (mod p)`, 즉 **N의 미지 약수 p(≈N^0.5)를 법으로 하는 다변수 선형식의 작은 근** 문제.
→ Herrmann–May, "Solving Linear Equations Modulo Divisors: On Factoring Given Any Bits" (ASIACRYPT 2008).

파라미터:
- β = 0.5 (p ≥ N^β).
- 변수 bound: X1=2^84, X2=2^58, X3=2^69, X4=2^87, X5=2^46.
- multiplicity m, extra shift t 는 튜닝. n=5에서 344비트를 잡으려면 m을 4~7 정도로 키워야 함(격자 차원 ↑).

---

## 5. 구현 방법 A — SageMath (권장)

Sage의 `small_roots`는 1변수 전용이라 다변수 divisor 케이스는 외부 구현을 쓴다.
가장 흔한 선택: **defund/coppersmith** (https://github.com/defund/coppersmith) 의 `coppersmith.sage`.
그 안의 `small_roots(f, bounds, m, d)` 를 그대로 사용.

```python
# sage 파일: solve.sage
load("coppersmith.sage")   # defund/coppersmith

# (2절의 N,e,ct,MASK,P_AND_MASK 정수들을 여기에 정의)

offsets = [265, 362, 600, 682, 784]
lengths = [84, 58, 69, 87, 46]

PR = PolynomialRing(Zmod(N), 'x', 5)
xs = PR.gens()

def try_pair(b150, b920, m=5, d=3):
    c = P_AND_MASK + (b150 << 150) + (b920 << 920)
    f = c
    for xi, off in zip(xs, offsets):
        f += xi * (1 << off)
    bounds = [2**L for L in lengths]
    # divisor 버전: beta=0.5 를 지원하는 구현을 쓰거나,
    # f 를 mod N 다변수로 두고 small_roots(f, bounds, m, d) 호출.
    roots = small_roots(f, bounds, m=m, d=d)
    for r in roots:
        p = c + sum(int(r[i]) * (1 << offsets[i]) for i in range(5))
        if p > 1 and N % p == 0:
            return int(p)
    return None

def run():
    for b150 in range(16):
        for b920 in range(16):
            p = try_pair(b150, b920)
            if p:
                print("FOUND p =", hex(p))
                return p
    print("no solution — increase m/d")

run()
```

> 주의: defund 구현은 "mod N" 다변수 Coppersmith다. divisor(β<1) 케이스를 제대로 쓰려면
> Herrmann–May 전용 시프트(예: `f^k · N^(max(t-k,0))`, `t = floor(m·β)` 근처)를 쓰는 구현이 필요.
> β=0.5 divisor를 지원하는 스크립트(예: "factoring with known bits" HM 구현)를 쓰는 것이 안전하다.
> 첫 시도는 작은 m(=4,5)로, 실패하면 m,t를 키운다. 격자 차원 ≈ C(m+5,5).

병렬화: `try_pair`는 256쌍이 독립 → `multiprocessing.Pool`로 분산, 정답 나오면 조기 종료.

---

## 6. 구현 방법 B — Sage 없이 (fpylll + sympy)

Sage가 없으면 직접 구성. 큰 틀:

1. `pip install fpylll sympy` (또는 `pip install --break-system-packages ...`).
2. **시프트 다항식 생성.** monomial 집합 `{ x1^a1 ... x5^a5 : Σai ≤ m }` 각각에 대해
   `g = (mono) · f^k · N^(max(t-k,0))` 형태로 divisor 시프트 구성 (t = round(m·β)).
   스케일링: 변수 xi → Xi·xi 치환해 계수 크기 정규화.
3. **격자 B** = 각 g의 (정규화된) 계수 벡터를 행으로. fpylll `LLL.reduction(B)`.
4. 축약된 짧은 행벡터 → Xi로 역스케일해서 정수계수 다항식 hi(x1..x5) 복원.
   이들은 참근을 정수 위에서 공유.
5. **연립 풀이.** sympy `groebner([h1,h2,...], x1..x5, order='lex')` 또는 순차 resultant로
   x1..x5 정수해 추출.
6. `p = c + Σ xi·2^offset_i`, `N % p == 0` 확인.

핵심 함수 시그니처 예:

```python
def herrmann_may(c, offsets, lengths, N, beta=0.5, m=5, t=None):
    # returns list of candidate (x1..x5) integer tuples
    ...
```

LLL은 fpylll로 충분(차원 수백까지 실용). 근 추출이 까다로우면 m을 키워 짧은 벡터를 더 확보.

---

## 7. p 획득 후 복호화 (표준 RSA)

```python
q = N // p
assert p * q == N
phi = (p - 1) * (q - 1)
d = pow(e, -1, phi)              # py3.8+
m = pow(ct, d, N)
pt = m.to_bytes((m.bit_length() + 7) // 8, "big")
print(pt)                        # 최종 평문 바이트열 (big-endian) = 정답
```

검증 포인트:
- `p*q == N` 필수. (원하면 p,q 소수성도 확인)
- 복호 결과 바이트열이 의미 있는(ASCII/플래그 형태) 값인지 확인.
- 평문이 짧으면 앞쪽 0바이트 패딩 주의(big-endian 그대로 출력).

---

## 8. 브루트포스 비용 메모

- 조합 수 자체는 256 → 무시 가능.
- 실비용 = 256 × (격자 1회 축약). 정답 시 조기 종료 → 평균 ≈128회.
- n=5, m 튜닝에 따라 1회 수십 초~수 분. 단일코어 순차 ≈ 수 시간, 16코어 병렬 ≈ 수 분~30분.
- 최적화: (1) m·t 최소화로 격자 축소, (2) 256쌍 멀티프로세스 병렬, (3) `N%p==0` 즉시 중단.

---

## 9. 실행 순서 체크리스트

1. 2절 데이터 로드 → 3절 파싱으로 blocks 재확인 (offset/length가 표와 일치하는지).
2. 격자 구현 선택 (A: Sage/defund, B: fpylll 직접).
3. 작은 m(=4)으로 1쌍(b150=b920=0) 테스트 → 파이프라인 동작 확인.
4. m을 올려가며 344비트가 잡히는 최소 m 확보.
5. 256쌍 루프(병렬) 실행 → p 획득.
6. 7절로 복호화 → 평문 바이트열 출력.
