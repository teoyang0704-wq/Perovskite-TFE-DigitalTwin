# TFE Physics-Informed Digital Twin — 문헌 조사 및 Research Gap 분석 (Step 0)

- 작성일: 2026-07-06 (웹 검색 기반 초안)
- 프로젝트 목표(확정 문안): **"기존 문헌에 있는 물성 데이터와 물리 모델을 활용하여, 페로브스카이트 태양전지용 TFE 설계를 체계적으로 탐색하고 설계 지침을 도출하는 Physics-Informed Digital Twin Framework를 개발한다."**
- ⚠️ 주의: 아래 인용은 검색 결과(초록·본문 발췌) 기반으로 정리한 것입니다. **논문 본문에 인용하기 전 반드시 원문(DOI)을 직접 확인**하세요. 각 항목에 확인 상태를 표기했습니다.
  - ✅ 서지사항(저자·저널·연도·DOI)까지 확인됨
  - 🔶 내용은 확인, 서지 일부(권·호·연도 등) 재확인 필요

---

## 0. 요약 (Executive Summary)

4개 축(A: 봉지 실험·리뷰 / B: 투습 물리 모델링 / C: 디지털 트윈 / D: AI 최적화·XAI)의 문헌을 교차 분석한 결과:

1. **봉지 실험**은 우수한 소재·공정 조합(예: ALD Al₂O₃ / iCVD 유기층 dyad)을 보고하지만, 기하 설계공간(두께 × 쌍수 × 순서)의 탐색은 소수 조건의 시행착오에 머문다.
2. **투습 물리 모델**(라미네이트 직렬저항, lag-time, 핀홀/굴곡경로)은 확립되어 있으나, 등온·정상상태 해석 도구에 머물며 동적 환경(온·습도 일주기), 열–수분 연성, 열기계 내구성, 수명 예측과 연결되지 않는다.
3. **PV 디지털 트윈**은 대부분 시스템/인버터/발전량 수준의 데이터 기반 모델이며, 봉지 스택(소자 레벨) 설계 탐색용 물리 기반 트윈은 부재하다.
4. **AI 최적화·XAI**는 (i) 광·전기 성능 단일 도메인의 두께 최적화이거나 (ii) 조성·재료 선택의 데이터 기반 ML+SHAP이며, 봉지의 다목적(수명·열·광학·비용·무게·내구성) 트레이드오프와 물리 인과 해석을 결합한 사례가 없다.

→ **Gap의 본질**: 개별 부품(물성 데이터, 1D 투습 모델, NSGA-II, SHAP)은 모두 존재하지만, 이들을 **"봉지 기하 설계"라는 문제에 대해 검증 가능하게 통합한 프레임워크**가 없다. 본 프로젝트의 기여는 바로 이 통합과 그로부터 나오는 정량 설계 지침이다.

---

## 1. 조사 범위와 방법

| 축 | 검색 키워드(대표) | 목적 |
|---|---|---|
| A. 봉지 실험·리뷰 | perovskite encapsulation, TFE, ALD Al2O3, WVTR, dyad | 실험 현황·정량 anchor 확보 |
| B. 투습 물리 모델 | multilayer barrier permeation, lag time, pinhole, tortuous path | 채택할 물리 모델 계보 확인 |
| C. 디지털 트윈 | digital twin photovoltaic, degradation prediction | DT 용어의 기존 용법과의 거리 확인 |
| D. AI·최적화·XAI | NSGA-II solar cell, machine learning perovskite stability, SHAP | 방법론 선행연구와 차별점 확인 |

---

## 2. 문헌 비교표

### A. 페로브스카이트/TFE 봉지 — 실험 및 리뷰

| 문헌 | 목적·방법 | 핵심 결과 | 한계 | 본 프로젝트와의 차이 |
|---|---|---|---|---|
| Lee et al., *Adv. Energy Mater.* (2017/18) — iCVD/ALD 다층 TFE의 PSC 적용 🔶 | pV3D3(iCVD)/Al₂O₃(ALD) 교대 dyad TFE를 PSC에 최초 적용; 공정온도·퍼지 최적화 | 3-dyad에서 WVTR ~10⁻⁴ g/m²/day (38 °C/90 %RH); **1–4 dyad의 WVTR·lag time 실측 곡선 제공** | dyad 수 1–4, 두께 소수 조건만 실험; 설계공간 지도·트레이드오프(광학/비용/무게/열) 미제시 | 우리는 동일 재료계의 (d_org, d_inorg, n) 연속 설계공간을 시뮬레이션으로 전수 탐색하고 6목적 파레토를 제시. **이 논문의 dyad–WVTR·lag 곡선은 Step 5 검증 타깃 V2** |
| Choi et al., *Sol. Energy Mater. Sol. Cells* 179, 80 (2018) 🔶 | 저온(95 °C) ALD Al₂O₃ 단층 봉지; HTM별 열화 분석 | 50 nm Al₂O₃ WVTR 1.84×10⁻² g/m²/day (45 °C/100 %RH); 봉지 소자 7500 h(RT/50 %RH) 후 PCE 감소 <4 % | 단층·단일 두께 중심; 다층 기하 효과, 열·기계·광학 목적 미고려 | 단층 두께–WVTR 관계는 **핀홀 폐색 모델 캘리브레이션(검증 타깃 V1)** 으로 사용 |
| Ramos et al., *Sustain. Energy Fuels* 2, 2468 (2018) ✅(제목·저널) | 60 °C 저온 ALD-Al₂O₃ 단일 코팅 봉지의 소자 호환성 | 저온 공정으로 PSC 손상 없이 장기 안정성 개선 | 위와 동일: 단층, 기하 변수 미탐색 | 공정 온도 제약을 DB의 '공정 제약 필드'로 반영 |
| Singh et al., *Sol. Energy Mater. Sol. Cells* 205, 110289 (2020) 🔶 | PSC 위 ALD Al₂O₃ 성장 메커니즘 + 장기 안정성 | 300일 간헐 측정에서 초기 효율 84 % 유지 | 성장 메커니즘 중심; 다층·다목적 설계 없음 | 계면 성장/결함 논의를 결함 모델 파라미터(f₀, d_close) 근거로 인용 |
| Wu(?) et al., OLED 다층 배리어 (PMC9078233, 2022) 🔶 | ALD Al₂O₃ + parylene C dyad OLED 봉지; 결함 원인 규명 | 단층 고품질 무기막(20–30 nm) 한계 ~10⁻⁴; **3 dyad(30 nm/500 nm)로 <10⁻⁵ g/m²/day; ALD 중 파티클이 핵심 결함원** | OLED 대상; 경험적 최적화; 모델 없음 | '파티클 잔류 결함 바닥(f_res)'의 실험 근거; **검증 타깃 V3** |
| Dipta, Rahim & Uddin, *Appl. Phys. Rev.* 11, 021301 (2024) ✅ | 봉지+납 누출 방지 리뷰 | 수분·Pb²⁺ 동시 차단 봉지 필요; **표준화된 안정성 시험 프로토콜 부재 지적** | 리뷰: 정량 설계 도구 부재를 스스로 확인 | 우리 결과 보고 시 시험조건(85/85, ISOS 등) 명시 근거 |
| *ACS Materials Au* 봉지·안정성 시험 리뷰 (2022) 🔶 | 봉지 구조·시험법 리뷰 | 안정성 리뷰 대비 봉지 리뷰 희소; 우수 봉지로 damp heat >800 h·옥외 1700 h 사례; ~200 h 가속 스크리닝 제안 | 실험 리뷰; 시뮬레이션 기반 설계 탐색 부재 | 서론에서 "봉지 연구의 병목 = 느린 실험 반복" 논거로 인용 |

### B. 투습 물리 모델링·시뮬레이션

| 문헌 | 목적·방법 | 핵심 결과 | 한계 | 본 프로젝트와의 차이 |
|---|---|---|---|---|
| **Graff, Williford & Burrows, *J. Appl. Phys.* 96(4), 1840–1849 (2004), DOI 10.1063/1.1768610** ✅ | 다층 유/무기 배리어의 투습 메커니즘을 과도·정상상태 측정 + Fick 확산 모델로 규명 | 다층 배리어 성능은 평형 확산이 아니라 **긴 유효 경로(굴곡 경로)에 의한 lag-time이 지배**; t_lag = l²/6D | 해석적 정상상태/lag 프레임; 동적 환경·열 연성·수명·다목적 없음 | 본 엔진의 굴곡경로(τ²)·핀홀 유효매질의 이론적 근원. 우리는 이를 **시간영역 1D PDE + Arrhenius + 동적 BC**로 확장 |
| Kiese et al., *Thin Solid Films* (2019) — 시간의존 다층 투습: 실험 vs 이론 🔶 | 이상 라미네이트(직렬저항) 모델 + 해석적 lag-time 식과 실측 비교 | 나노결함 준균질 투과 가정 하 직렬저항 모델 유효; B/P 적층 조합론적 lag 식 | 등온; 설계 탐색·최적화 없음 | 우리 트윈의 **정상상태 해는 동일 직렬저항 해석해와 기계 정밀도로 일치(자체 검증 완료)**; 그 위에 과도·동적 해석을 얹음 |
| 무기막 핀홀·입계 3D 확산 시뮬레이션 계열 (PET 위 단/다층; 서지 재확인 필요) 🔶 | 핀홀·입계를 명시적으로 넣은 3D 수치 확산 | WVTR가 무기층 두께에 지수적으로 감소; 결함 간격·입계가 지배 | 계산 비용 큼 → 최적화 루프 삽입 곤란; 다목적 없음 | 우리는 3D 결함 효과를 **1D 유효매질(f_pin, τ²)로 축약**해 1회 평가 ~1 s를 달성, 수천 회 평가가 필요한 NSGA-II에 삽입 가능 |

### C. 디지털 트윈 (PV 분야에서의 용법)

| 문헌 | 목적·방법 | 핵심 결과 | 한계 | 본 프로젝트와의 차이 |
|---|---|---|---|---|
| Zhang et al., *J. Mod. Power Syst. Clean Energy* 12(5), 1472 (2024) 등 PV 발전량 DT ✅(서지 문자열) | 실시간 데이터 동기화 기반 발전량 예측 DT | 데이터 기반 예측·이상 감지 | **시스템 레벨**; 소자 물리·설계 변수 없음 | 우리는 '설계 탐색용 소자(봉지 스택) 레벨 물리 트윈' — DT 용어의 다른 용법임을 서론에서 명확히 구분 |
| 인버터 열화 추정 ML-DT (IEEE ECCE 2024, DOI 10.1109/ECCE55643.2024.10861035) ✅(서지 문자열) | 회로 레벨 DT + ML로 인버터 부품 열화 추정 | 부품 파라미터 열화 추적 | 전력전자 부품 대상 | 동일: 레벨과 목적이 다름 |
| Gr–Si 셀 physics-informed semi-empirical DT, *Sci. Rep.* (2026) 🔶 | 방사선 열화 데이터(합성 포함) + Random Forest로 소자 열화 예측 | R² > 0.96 예측 | 회귀 대리모델 중심; 봉지·수분·설계 탐색 아님; 합성데이터 검증 | 소자 레벨 'physics-informed DT' 용어의 선례로 인용하되, 우리는 **PDE 순방향 모델 자체가 트윈**이고 문헌 실측으로 검증한다는 점을 대비 |

### D. AI 기반 최적화·해석 (ML / MOO / XAI)

| 문헌 | 목적·방법 | 핵심 결과 | 한계 | 본 프로젝트와의 차이 |
|---|---|---|---|---|
| **Hartono et al., *Nat. Commun.* 11, 4172 (2020), DOI 10.1038/s41467-020-17945-4** ✅ | 21종 유기할라이드 캡핑층을 특징화 → 가속 열화 실험 → RF 회귀 + SHAP로 안정성 지배 인자 도출, 재료 설계 지침 생성 | 수소결합 공여자 수·극성표면적이 안정성과 상관; PTEAI가 수명 4±2배 연장 | **조성/재료 선택** 문제; 데이터 기반(물리 모델 없음); SHAP는 상관 수준 해석; 기하 변수 없음 | 우리는 **기하 변수 + 물리 PDE 순방향 모델** 기반. 해석도 SHAP 수치 → 핀홀 폐색/굴곡경로/균열이라는 **검증 가능한 인과 사슬**로 번역(우리 Step 7–8). 가장 직접적인 방법론적 대조군 |
| GaAs 박막 셀 NSGA-II 최적화, *Solar Energy* (2021) 🔶 | 활성층+3층 ARC 두께를 NSGA-II로 흡수↑·비용↓ 파레토 탐색 | 흡수 50 %·Jsc 43 %·효율 45 % 개선점 도출 | **광·전기 성능 단일 도메인**; 신뢰성(수분·열·기계) 목적 없음 | 우리는 신뢰성 중심 6목적(수명·열·광투과·비용·무게·내구성)으로 확장; 층두께 MOO의 방법적 선례로 인용 |
| Moulebhar et al., *Phys. Status Solidi A* (2025), DOI 10.1002/pssa.202400654 ✅ | RSM(BBD)+NSGA-II+SCAPS-1D로 유기 셀 층두께·도핑 최적화 | 최적 두께/도핑 조합 제시 | 효율 중심; 봉지·수명 없음 | 동일 취지의 대조군: '시뮬레이터+NSGA-II'는 확립된 조합 → 우리의 novelty는 알고리즘이 아니라 **문제 정의(봉지 기하×6목적)와 검증·해석 체계**임을 명시 |
| ARC 다층 반사율 DL 예측, *Opt. Quantum Electron.* (2025) 🔶 | TMM 데이터로 학습한 DL 대리모델 | 고속 반사율 예측 | 광학 단일 물성 | 대리모델 접근의 참고; 우리는 물리 엔진이 충분히 빨라(≈1 s/평가) 대리모델 없이 직접 최적화 가능 |

---

## 3. Research Gap 종합

- **G1. 기하 설계공간의 체계적 탐색 부재** — 실험(A군)은 dyad 수 1–4, 두께 몇 조건 수준의 비교에 그친다. (두께_유기 × 두께_무기 × 쌍수 × 순서)의 연속 설계공간을 정량 지도로 그린 연구가 없다.
- **G2. 다물리·동적 환경 통합 부재** — B군 모델은 등온 정상상태/lag 해석에 강하지만, 온·습도 일주기 + Arrhenius D(T) + 열전도 + 열기계 응력 + 누적 수분량 기반 수명(T80)을 하나의 시간영역 시뮬레이션으로 잇는 프레임워크가 없다.
- **G3. 목적함수의 단일성** — D군 MOO는 광·전기 성능/비용에 국한. 봉지가 실제로 직면하는 6대 상충 지표를 파레토 구조로 다룬 연구가 없다.
- **G4. 해석의 통계성** — SHAP 기여도(Hartono 2020)는 "무엇이 중요했는가"까지만 답한다. "왜 그런가"를 물리 법칙(핀홀 폐색, 굴곡 경로, 임계 균열 두께)과 매핑해 설계 규칙으로 자동 번역하는 체계가 없다.
- **G5. 검증 절차의 부재** — DT/ML 연구 다수가 자체·합성 데이터로 검증한다. 문헌의 독립 정량 데이터(단층 두께–WVTR, dyad 수–WVTR·lag)에 대한 교차 검증을 명시적 단계(우리 Step 5)로 두는 연구가 드물다.

> 요약: **부품은 다 있으나 조립된 적이 없다.** 조립 자체와, 조립이 처음으로 가능하게 하는 산출물(6목적 파레토 지도 + 물리 인과 기반 설계 규칙)이 본 연구의 자리다.

---

## 4. Research Question

**Main RQ.** 문헌 보고 물성과 확립된 결함 매개 1D 다물리 모델로 구성·검증한 디지털 트윈은, 페로브스카이트 TFE의 기하 설계공간에서 6대 성능 지표의 트레이드오프 구조와 정량적 설계 규칙을 실험 문헌과 일관되게 도출할 수 있는가?

- **RQ1 (충실도/검증).** 캘리브레이션된 트윈이 (i) 단층 무기막 두께–WVTR 경향과 (ii) dyad 수(1–4)–WVTR·lag time 경향을 실험 산포(±반 자릿수 수준) 내에서 재현하는가?
- **RQ2 (트레이드오프 구조).** 6목적 파레토 프론트에서 기하 변수별 지배 상충관계는 무엇이며(예: 쌍수↑ → T80↑ vs 투과율·비용·무게 저하), knee-point 설계는 어디에 형성되는가?
- **RQ3 (설계 규칙과 인과).** 도출된 규칙(무기층 성능 포화 두께, 균열 임계 이하 유지, 최적 쌍수 범위, 광투과 유지를 위한 총두께 상한)을 핀홀 폐색–굴곡 경로–열기계 균열의 물리 인과 사슬로 설명하고, 문헌 관찰과 대조해 정당화할 수 있는가?

---

## 5. Novelty (겸손한 프레이밍)

**주장하지 않는 것(Not-claims):** 새로운 재료 ✗ / 새로운 PDE·확산 이론 ✗ / 새로운 최적화 알고리즘 ✗ / "세계 최초 디지털 트윈" ✗

**주장하는 것(Claims): 통합과 체계화, 그리고 그 산출물**
1. **문제 재정의**: 봉지 연구를 '조성 탐색'이 아닌 **기하 설계 문제**로 정식화 (설계변수: d_org, d_inorg, n_pairs, 적층 순서).
2. **경량 결함 매개 1D 다물리 트윈**: Graff류 굴곡경로·핀홀 물리를 1D 유효매질로 축약해 시간영역 PDE(수분+열, Arrhenius 연성, 동적 온습도 경계조건)로 확장 — 1회 평가 ~1 s로 최적화 루프 삽입 가능(상용 3D FEM 대비 실용적 차별점).
3. **6목적 파레토 + 물리 매핑 해석**: 신뢰성 중심 6목적 동시 최적화와, SHAP 수치를 물리 인과 사슬로 번역하는 설계 지침 자동 생성.
4. **문헌 교차 검증 절차의 명시화**: 독립 실험 정량 데이터(아래 §8 V1–V4)를 검증 타깃으로 고정.

한 줄 포지셔닝: *"우리는 새로운 물리를 만들지 않는다. 흩어져 있던 검증된 물리·데이터·알고리즘을 봉지 기하 설계 문제 위에 처음으로 조립하고, 그 조립이 만들어내는 설계 지도를 제시한다."*

---

## 6. Expected Contribution

1. **(방법론)** 재사용 가능한 오픈 프레임워크: 문헌 물성 DB 스키마 + 경량 1D 결함매개 PDE 엔진 + NSGA-II 6목적 최적화 + 물리 매핑 리포트 생성기 (코드·DB 공개 시 커뮤니티 기여).
2. **(과학적 이해)** 기하 변수 → 6대 지표의 정량 지도와 지배 메커니즘 규명: 실험이 산발적으로 보고한 경향(dyad 증가 시 WVTR·lag 개선, 두꺼운 단층의 균열 취약성)을 하나의 인과 프레임으로 통합 설명.
3. **(실용 설계 지침)** 공정 개발자가 즉시 쓸 수 있는 규칙 형태의 산출물: 예) "무기층 X nm 이상에서 핀홀 폐색 포화", "쌍수 N₁–N₂가 수명–광학–비용 knee", "총두께 Y nm 이하에서 투과율 Z 이상 유지" (수치는 시뮬레이션+검증 후 확정).
4. **(검증 문화)** 시뮬레이션 기반 설계 연구에서 문헌 정량치 교차검증을 표준 절차로 제시하는 템플릿.

---

## 7. 예상 리뷰어 공격과 방어 (사전 점검)

| 예상 지적 | 방어 전략 |
|---|---|
| "1D로 3D 결함 물리를 담을 수 있나?" | 결함 매개 유효매질(f_pin, τ²)의 이론적 근거(Graff 2004 계보) 제시 + Step 5에서 실험 경향 재현으로 실증 + 한계(측면 불균일, 국소 결함 통계) 명시 |
| "물성값 불확실성" | Step 1 DB에 값의 범위·측정조건·출처를 필드로 기록; 민감도 분석(모든 KPI에 대해) 수행 |
| "비용·무게 모델이 단순" | 프록시(공정시간·재료량 비례)임을 명시하고 상대 비교 목적으로 한정 |
| "T80 임계 수분량 M_crit의 자의성" | 캘리브레이션 파라미터로 명시, 문헌 damp-heat 수명과 역산 대조, 민감도 제시 |
| "DT라는 용어 남용 아닌가" | §2C처럼 기존 용법(시스템 모니터링 DT)과 구분: '설계 탐색용 오프라인 물리 트윈'으로 정의 |

---

## 8. 10단계 로드맵 ↔ 문헌 연결

### Step 1 (물성 DB) 필드 후보 ← A·B군
Material(계열/공정) · D(T_ref), Ea 또는 WVTR@(T,RH) · 용해도 S · k(열전도) · ρ · c_p · n(굴절률) · α_opt · CTE · E(영률) · ν · 임계 균열 두께 · 두께 범위 · 공정 온도 제약 · 비용 프록시 · 출처(DOI)·측정 조건 — *구체 스키마는 다음 작업(Step 1)에서 확정.*

### Step 5 검증 타깃 (정량 anchor)
| ID | 데이터 | 출처 | 용도 |
|---|---|---|---|
| V1 | 단층 ALD Al₂O₃ 두께–WVTR (예: 50 nm → 1.84×10⁻² @45 °C/100 %RH; 고품질 20–30 nm → ~10⁻⁴) | Choi 2018 🔶 / OLED 2022 🔶 | 핀홀 폐색 파라미터(f₀, d_close, f_res) 캘리브레이션 |
| V2 | dyad 수 1–4 vs WVTR·lag time 곡선 (3-dyad ~10⁻⁴ @38 °C/90 %RH) | Lee 2017/18 🔶 | 다층 굴곡경로(τ²)·쌍수 효과 검증 (핵심) |
| V3 | 3 dyad (30 nm Al₂O₃/500 nm parylene C) < 10⁻⁵ g/m²/day; 파티클 결함 지배 | OLED 2022 🔶 | f_res(파티클 바닥) 근거·이종 재료계 일반화 확인 |
| V4 | 다층 배리어의 lag-time 지배 체제; PET 다층 <10⁻⁵ @25 °C/40 %RH | Graff 2004 ✅ | 과도 응답(파과 곡선) 정성·정량 검증 |
| V5 | 시간의존 투습 실험–이론 비교 데이터 | Kiese 2019 🔶 | 직렬저항+lag 해석해 대비 수치해 검증 보강 |

### Step 6–8 방법 선례
NSGA-II 층두께 최적화(Solar Energy 2021; pss(a) 2025) · SHAP 기반 설계 지침(Hartono 2020) → 알고리즘은 차용, **목적함수 구성과 해석 체계가 차별점**임을 논문에 명시.

---

## 9. 참고문헌 (확인 상태 포함)

1. ✅ G. L. Graff, R. E. Williford, P. E. Burrows, "Mechanisms of vapor permeation through multilayer barrier films: Lag time versus equilibrium permeation," *J. Appl. Phys.* 96(4), 1840–1849 (2004). DOI: 10.1063/1.1768610
2. ✅ N. T. P. Hartono et al., "How machine learning can help select capping layers to suppress perovskite degradation," *Nat. Commun.* 11, 4172 (2020). DOI: 10.1038/s41467-020-17945-4
3. ✅ S. S. Dipta, M. A. Rahim, A. Uddin, "Encapsulating perovskite solar cells for long-term stability and prevention of lead toxicity," *Appl. Phys. Rev.* 11, 021301 (2024). DOI: 10.1063/5.0197154
4. 🔶 Lee et al., "A Low-Temperature Thin-Film Encapsulation for Enhanced Stability of a Highly Efficient Perovskite Solar Cell," *Adv. Energy Mater.* (2017/2018) — 권·호·페이지 확인 필요
5. 🔶 Choi et al., "Enhancing stability for organic-inorganic perovskite solar cells by atomic layer deposited Al₂O₃ encapsulation," *Sol. Energy Mater. Sol. Cells* (2018) — S0927024818304239
6. ✅ Ramos et al., "Versatile perovskite solar cell encapsulation by low-temperature ALD-Al₂O₃ with long-term stability improvement," *Sustain. Energy Fuels* 2 (2018)
7. 🔶 Singh et al., "ALD Al₂O₃ on hybrid perovskite solar cells: Unveiling the growth mechanism and long-term stability," *Sol. Energy Mater. Sol. Cells* (2020) — S092702481930618X
8. 🔶 Kiese et al., "Time-dependent water vapor permeation through multilayer barrier films: Empirical versus theoretical results," *Thin Solid Films* (2019) — S004060901930001X
9. 🔶 "Efficient multi-barrier thin film encapsulation of OLED using alternating Al₂O₃ and polymer layers" (2022) — PMC9078233, 저널·저자 확인 필요
10. ✅(서지 문자열) X. Zhang et al., "Digital twin empowered PV power prediction," *J. Mod. Power Syst. Clean Energy* 12(5), 1472–1483 (2024)
11. ✅(서지 문자열) S. Roy et al., "An ML-enhanced digital twin model of photovoltaic inverter for estimating component degradation," IEEE ECCE 2024. DOI: 10.1109/ECCE55643.2024.10861035
12. 🔶 "Digital twin analysis of graphene–silicon solar Schottky-heterojunction cell efficiency," *Sci. Rep.* (2026)
13. 🔶 "Boosting photoelectric performance of thin film GaAs solar cell based on multi-objective optimization," *Solar Energy* (2021) — S0038092X21009907
14. ✅ Moulebhar et al., "Hybrid Optimization Approach Using Multiobjective Genetic Algorithm NSGA-II, SCAPS-1D Simulation, and Response Surface Methodology for Organic Solar Cell Analysis," *Phys. Status Solidi A* (2025). DOI: 10.1002/pssa.202400654
15. 🔶 "Encapsulation and Stability Testing of Perovskite Solar Cells for Real Life Applications," *ACS Materials Au* (2022). DOI: 10.1021/acsmaterialsau.1c00045
16. 🔶 무기막 핀홀·입계 3D 확산 시뮬레이션 (PET 기판) — 서지 전체 확인 필요

*다음 단계(Step 1): 위 A·B군 문헌에서 실제 물성값을 수집할 CSV/DB 스키마 설계 — 코드 없이 데이터 구조 정의부터.*
