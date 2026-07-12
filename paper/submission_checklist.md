# 투고 마스터 체크리스트 — 사람이 해야 할 일 전부 (단일 문서)

작성일 2026-07-07 · 이 문서 하나만 보고 작업하면 됩니다. 완료 시 [ ] → [x].

---

## A. 서지 확인 ([VERIFY] 12건) — 원문 DOI로 대조 후 원고 v2 참고문헌 갱신
- [x] ref 2 확정: Wang et al., ACS Mater. Au 2 (2022) 215–236 ~~ACS Mater. Au 2022 리뷰 — 저자·권·페이지
- [x] ref 3 확정: Y. I. Lee et al., AEM 8 (2018) 1701928 (+외부 프로브 SI S4.5 완료) ~~Lee — 권·호·페이지 (★Fig.2d 데이터 논문)
- [x] Wu(7)·Kempe(10)·Dameron(20) 원문/참고문헌으로 재검증 완료(2026-07-10); Meyer 2009 AM 21,1845 / Affinito 2000 SCT 133,528 / Moro 2004 SPIE 5214,83 추가 가용 인용 확보
- [x] ref 4 확정: Choi·Kim·Lim·Han·Ho-Baillie·Park, SEMSC 188 (2018) 37–43 (페이지만 1분 재확인 권장; RH 모순 SI 반영) ~~Choi — 전체 서지 + 측정법(§BM013 quality B→A)
- [x] ref 5 확정: Ramos et al., Sustain. Energy Fuels 2 (2018) 2468–2479 ~~Ramos et al., Sustain. Energy Fuels 2018 — 페이지
- [x] ref 6 확정: Singh·Ghosh·Subbiah·Mahuli·Sarkar, SEMSC 205, 110289 (2020) ~~ref 6 Singh et al., SEMSC 205, 110289 — 확인만
- [x] ref 9 확정: Kiese·Kücükpinar·Miesbauer·Langowski, TSF 672 (2019) 199–205 ~~ref 9 et al., Thin Solid Films 2019 — 저자·권·페이지
- [x] ref 12 확정: Zhang·Qi·Ji·He·Ren·Li, Sol. Energy 230 (2021) 1122–1132 ~~ref 12 NSGA-II, Solar Energy 2021 — 저자·권
- [~] ref 14: 제목·저널 확정(JCP, PII 기재; 초록이 §1 차별화 뒷받침 확인) — 저자·권만 1분 ~~ref 14 PSC DT, J. Clean. Prod. 2024 — 저자·권 (★정독 대상)
- [~] ref 15: 제목·저널 확정(Mater. Sci. Eng. B 2025, PII 기재) — 저자·권만 1분 ~~ref 15 2025 — 저널 확정(Thin Solid Films?)·저자 (★정독 대상)
- [x] ref 17 확정: JAP 74(9) 5471–5475 ~~ref 17 & Nulman, J. Appl. Phys. 74 (1993) — 권·페이지
- [x] ref 18 확정: JVST A 18(1) 149–157 ~~ref 18 Sobrinho, JVST A 18 (2000) — 페이지
- [x] ref 22 확정: JPC 61(1) 93–95 ~~ref 22, J. Phys. Chem. 61 (1957) — 페이지

## B. 근접 논문 정독 (★4편) — 각각 "1문장 차별"이 원고 §1과 일치하는지 확인
- [ ] AutoML PSC digital twin (2024, S0959652624040095)
- [ ] TFE figure-of-merit (2025, S0921510725008451)
- [ ] Hartono et al. 2020 (서지는 OK, 내용 정독)
- [ ] Dual-model ML efficiency+stability (Chem. Eng. J. 2025)
- [ ] (+ Graff 2004 원문에서 lag/tortuous 논지 및 vertical–horizontal 등가 확인 — SI S3 근거)

## C. 과학 잔여 항목
- [ ] **72 h 출처 재확인 필요**: 본문 전문에 시험시간 언급 없음(2026-07-10 확인). 어디서 보셨는지 확인 — MOCON 매뉴얼? 타 논문? 원고는 지속시간-비의존(2.0–5.5일 창)으로 전환 완료라 투고엔 지장 없음
- [x] **M_crit 캘리브레이션 완료(2026-07-07)**: Tier-1 하한 0.97 g/m² [0.90–1.04] (SI S4.4b) → 창 내 전 설계 수분수명 ≥25 yr, 수분은 더 이상 지배 채널 아님; 양측 추정은 향후 임계 초과 노화 데이터 필요
- [ ] **Lee Fig.2d 확보 시**: 디지타이저(data/provenance/digitize_fig3.py 재사용) → 독립 외부검증 실행 → §3.5(v) 완화. 저자 요청 이메일 초안 필요하면 요청
- [ ] 역학 1차문헌 2편 확보(Al₂O₃ 박막 E·σ_c·임계두께) → P0004/P0006 quality C 교체, G3 신뢰도 LOW-MED→MED
- [ ] (선택) parylene-C D·S 1차 문헌 → E012/E013 승격

## D. 저자·정책
- [ ] 저자명·소속·교신 이메일 확정 (원고·커버레터·레포 3곳)
- [ ] CRediT 문구 작성
- [ ] 대상 저널의 생성형 AI 정책 원문 확인 → Declarations 문구 최종화 (현재 템플릿: Claude 사용 고지 + 저자 책임)
- [ ] 이해상충·펀딩 문구

## E. 레포 공개 전
- [ ] 라이선스 선택(MIT 또는 Apache-2.0 권장) → LICENSE 파일 추가
- [ ] [repository URL] 실제 주소 삽입 — 원고 Declarations / SI S6 / 커버레터 3곳
- [ ] README에 인용 안내(CITATION.cff) 추가
- [ ] 업로드 전 최종 실행 테스트: S6의 6개 명령어를 클린 환경에서 1회

## F. 포맷 (마지막 단계 — 합의대로 최후에)
- [x] LaTeX 변환 완료: manuscript_v3.tex(+PDF 11pp, article 클래스·Times) — 투고 시 elsarticle로 1줄 교체
- [ ] Highlights 별도 파일 분리, 그림 개별 업로드(현재 300 dpi PNG + 벡터 PDF 준비됨), SI PDF화
- [ ] 커버레터의 [n]·추천 리뷰어·날짜 기입

## G. 산출물 지도 (현재 위치)
| 산출물 | 경로 |
|---|---|
| 원고 v2 (리뷰 반영) | manuscript/manuscript_draft_v2_SEMSC.md |
| 가상 리뷰 + 대응표 | manuscript/peer_review_and_response.md |
| 커버레터 | manuscript/cover_letter_SEMSC.md |
| SI (S1–S6) + Fig S1–S3 + Table S1 | SI/ |
| 본문 Fig.1–7 (PNG 300dpi + PDF) | paper_figures/ |
| 물성 DB (6 테이블 + provenance) | step1_db/ |
| 엔진·캘리브레이션·최적화·XAI·UQ 코드/리포트 | step5_calibration/ · step6_optimization/ · step7_xai/ · step8_uq/ |
| 문헌 리뷰·점유 분석 | TFE_literature_review_research_gap_v1.md · TFE_prior_art_novelty_analysis.md |
| GitHub형 재현 패키지(전체 사본) | repo_skeleton/Perovskite-TFE-DigitalTwin/ |

**투고 게이트 요약**: A(서지) + B(★정독) + C의 72h 위치·M_crit 두 개가 닫히면 F(포맷) 후 제출 가능. Lee 검증·역학 문헌은 "있으면 강함, 없어도 §3.5에 공개된 한계"로 처리되어 있음.

- [ ] **라이브 앱 데이터 갱신**: grid_all.csv·step6이 정정판(−6.7% 스케일)으로 재생성됨 — 다음 push에 포함 필요 (표시값 ~7% 차이, 임계값·순위 동일)

- [x] 서지 검증 세션 완료(2026-07-10): 22건 중 19건 완전확정, 3건(12·14·15)은 PII 앵커로 저자·권만 잔여(각 1분); 미인용 ref(pvdt2024) 제거

- [x] 외부검증 확장(2026-07-10): Carcia f-비교(0.09–0.11× floor)·Ea 교차(52.5 vs 40±10, +1.3σ)·Meyer 경계(×0.54)·Graff lag-기전 실험확인 — §3.1 한 문장 + SI S4.5(a–e) + CSV X005–X008
- [x] ~~Lee 데이터 요청~~ 불필요: 사용자가 원문 확보(2026-07-10) → **완전 모드-B 검증 완료** (−0.16/−0.34/−0.49 dec; SI S4.5a) — external_validation/author_request_email.txt (교신 S.G. Im, KAIST; 주소는 랩 페이지에서 확인 후 발송)

- [x] ★정독 Lee: 원문 전문 입수·정독 완료(구조·수치·방법 전부 반영) — 잔여 정독: FOM 2025, Hartono

- [x] **서지 22/22 완봉 (2026-07-10)**: ref 14 = Y. Yang, G. Jia, JCP 486 (2025) 144560 · ref 15 = M. Wu, M. Jian, L. Yang, MSE-B 323 (2026) 118821 — 사용자 스크린샷 + 인용 DB 교차 확정
- [x] AI 고지 정밀화(2026-07-10): 역할을 '구현·초안 보조'로 한정 + 소유권 문장 전면 배치 — 원고 Declarations·README·앱 푸터 일괄 정합(ISEF 팩 Q19/Q22와 일치)

- [x] elsarticle: 샌드박스에 클래스 없음(네트워크 차단) — **초회 투고는 Elsevier YPYW(Your Paper Your Way) 정책상 현 PDF 그대로 유효**; elsarticle 조판은 개정 단계에서만 필요(tex에 1줄 교체 주석 준비됨, PC에서 texlive-publishers 설치 후 컴파일)
