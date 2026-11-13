# AER-48 封存复核材料

这是完整的构造实验档案，用于复核方法，不声称来自真实厂站。48期各为6小时；on=1启用新曝气策略，0沿用原策略。主要问题是在登记模型下，能否确认长期每期耗电下降，而非比较两组原始均值。

## 观测材料
data/protocol.json是试验前登记的设计及参数格点；hot为试验前固定的高负荷标记。data/schedule.json是实际执行分配。data/meter.jsonl为按导出到达顺序保存的电表修订；period从1起，每期仅最高rev有效，同版冲突拒绝，低版迟到忽略。kwh是精确整数耗电；range=[L,U]表示所有闭区间整数仍可能，不能用已被高版撤销的旧读数填补。data/recovery.json是独立留存计数器对指定观测期的精确合计；所有记录同时成立。source.load与source.completions给出清洗后的材料和按完整y向量字典序编号的联合补全。报告须引用文件及JSON键或行号，区分观测、模型假设与推断。

## 实际随机化
每6期恰有3期on；同值连续不超过3期；全部hot期中恰有hot_on期on。每步先剔除无法延伸为完整合法分配的选择。只剩一个选择则概率1；两者均可时，若此前3期实际耗电和>=threshold，P(on)=2/3，否则1/3。最初缺少的历史用pre_y补齐，试验前on=0。此为在线分配，不是从合法分配中均匀抽样。design.Design保存调度源码，可用于核对单条路径；tiny_paths仅限不超过10期的小例，不能用于本批。

## 登记推断
参数tau、gamma分别从protocol中的两个格点数组取值，二者代表同期与滞后一期的整数kWh效应。登记的恒定加性模型为Y_i(a)=b_i+tau*a_i+gamma*a_(i-1)，a_0=0，b_i不随分配改变；这只是本次因果解释的前提，不得宣称已证明模型符合现实。对每份联合补全y及每对参数，由实际on反推出b，再让调度器在每条反事实路径上使用该路径产生的Y更新倾向。恢复合计只限制观测y的补全，不筛选反事实路径，也不能固定观测路径的概率后置换。

对任意合法反事实分配a，S(a)=sum(a_i*b_i)，其路径概率等于design.path_mass(a,b,tau,gamma)/3**n。令D=3**n，H(s)为得分s的整数路径质量和，M=sum(s*H(s))，s_obs=S(实际on)。精确双侧p值为sum(H(s) for abs(s*D-M)>=abs(s_obs*D-M))/D，等号计入尾部，禁止浮点阈值比较或蒙特卡洛替代。每格对全部联合补全取最大p；最大值>1/20才进入联合95%置信集合，同分见证取最小补全ID。该集合可不连通；对tau+gamma作精确集合投影，不能以包络区间填洞。只有投影非空且每个值<0才能在登记模型下确认持续节能；空集合报告材料与格点模型不相容，不作节能判定。不得把单个最优拟合或独立参数区间当联合结论。

## 本轮交付与入口
保留原始材料、source.py、design.py及checks.py原断言，只交付本批分析、REPORT.md、result.json和analysis.py。禁止枚举完整48期分配；须以合并等价前缀等方式累计全部路径的精确整数质量，说明何种历史信息足以同时决定后续可行性、倾向和检验得分，并论证没有错误合并。报告解释修订、补全、滞后与自适应分配如何改变结论，给出置信覆盖及最坏补全完整性的依据，记录实测耗时、峰值内存和适用边界。

从根目录调用python3 analysis.py run重算并覆盖REPORT.md和result.json；verify重新计算并核对报告及结构化证据，保留checks.py normal，并以不超过10期的小设计全枚举独立核对质量分布、尾部及补全极值，不能只读回结果判成功；demo打印本批实际联合集合、总效应支持集、极值见证与节能结论。更新.pointwise-runtime.json的verify/demo分别调用analysis.py verify/demo，setup仍为argv数组列表。

result.json为对象：cells按(tau,gamma)递增，每项含tau、gamma、p_by_completion（与source.completions顺序一致的约分分数[num,den]列表）、worst_id、p_max（同形分数）和accepted布尔；region为按(tau,gamma)递增的[参数,参数]列表；total_effect_support为去重递增整数列表；decision为saving、inconclusive或incompatible。附evidence对象，内容由分析者设计，须足以重跑核对递推、所有补全极值及报告数字，不以求解库成功状态代替证明。无合法补全时停止推断，报告相互冲突的材料位置，不输出部分置信集合。

当前运行清单只核对读取与调度的正常业务回归，尚未做任何效应推断；python3 checks.py demo展示实际修订与补全。最终verify/demo必须展示新增分析的真实结果。
