## 프로그램 설명
개인적으로 직접 개발한 CFD Solver 입니다.

외부 라이브러리는 intel MKL을 사용하고 있습니다.


## branch별 설명
### master
초기 파일들만 저장되어있습니다.

### template_version
solver에서 지원하는 모든 기능이 구현되어 있으며 template을 기반으로 코드가 작성되어있습니다.

### RDD
책임주도설계(responsibility driven design, RDD)에 대해 공부한 후 template_version branch의 코드를 객체지향에 맞게 리팩토링 중인 코드입니다.  
현재는 DG기반의 HOM에 대한 기능까지 리팩토링 되어 있습니다. (22.02.13)
