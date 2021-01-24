from pde import *
import warnings
warnings.filterwarnings('ignore')


class SGA:  # 最外层
    def __init__(self, num, depth, width, p_var, p_mute, p_rep, p_cro):
        # num: pool里PDE的数量
        # depth: 每个PDE的term的最大深度
        # width: 每个PDE所含term的最大数量
        # p_var: 生成树时节点为u/t/x而不是运算符的概率
        # p_rep: 将（所有）pde某一项重新生成以替换原项的概率
        # p_mute: PDE的树结构里每个节点的突变概率
        # p_cro: 不同PDE之间交换term的概率
        self.num = num
        self.p_mute = p_mute
        self.p_cro = p_cro
        self.p_rep = p_rep
        self.eqs = []
        self.mses = []
        self.ratio = 5
        print('Creating the original pdes in the pool ...')
        for i in range(num*self.ratio):
            # 循环产生num个pde
            a_pde = PDE(depth, width, p_var)
            a_err, a_w = evaluate_mse(a_pde)
            while a_err < 0.01 or a_err == np.inf:  # MSE太小则直接去除，to avoid u d t
                a_pde = PDE(depth, width, p_var)
                a_err, a_w = evaluate_mse(a_pde)
            # print('Creating the ith pde, i=', i)
            # print('a_pde.visualize():',a_pde.visualize())
            # print('evaluate_mse(a_pde):',a_err)
            self.eqs.append(a_pde)
            self.mses.append(a_err)

        new_eqs, new_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        sorted_indices = np.argsort(new_mse)
        for i, ix in enumerate(sorted_indices):
            self.mses[i], self.eqs[i] = new_mse[ix], new_eqs[ix]
        self.mses, self.eqs = self.mses[0:num], self.eqs[0:num]

    def run(self, gen=100):
        for i in range(1, gen+1):
            self.cross_over(self.p_cro)
            self.change(self.p_mute, self.p_rep)
            best_eq, best_mse = self.the_best()
            print('{} generation best_aic & best Eq: {}, {}'.format(i, best_mse, best_eq.visualize()))
            print('best concise Eq: {}'.format(best_eq.concise_visualize()))
            if best_mse < 0:
                print('We are close to the answer, pay attention')

    def the_best(self):
        argmin = np.argmin(self.mses)
        return self.eqs[argmin], self.mses[argmin]

    def cross_over(self, percentage=0.5):
        def cross_individual(pde1, pde2):
            new_pde1, new_pde2 = copy.deepcopy(pde1), copy.deepcopy(pde2)
            w1, w2 = len(pde1.elements), len(pde2.elements)

            ix1, ix2 = np.random.randint(w1), np.random.randint(w2)
            new_pde1.elements[ix1] = pde2.elements[ix2]
            new_pde2.elements[ix2] = pde1.elements[ix1]
            return new_pde1, new_pde2

        # 一般的好样本保存，并在此基础上交叉生成一半新样本
        # print('begin crossover')
        num_ix = int(self.num * percentage)
        new_eqs, new_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        sorted_indices = np.argsort(new_mse)
        for i, ix in enumerate(sorted_indices):
            self.mses[i], self.eqs[i] = new_mse[ix], new_eqs[ix]
        self.mses, self.eqs = self.mses[0:num_ix], self.eqs[0:num_ix]

        new_eqs, new_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        reo_eqs, reo_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        random.shuffle(reo_mse)
        random.shuffle(reo_eqs)

        for a, b in zip(new_eqs, reo_eqs):
            new_a, new_b = cross_individual(a, b)
            self.eqs.append(new_a)
            a_err, a_w = evaluate_mse(new_a)
            self.mses.append(a_err)

            self.eqs.append(new_b)
            b_err, b_w = evaluate_mse(new_b)
            self.mses.append(b_err)

        new_eqs, new_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        sorted_indices = np.argsort(new_mse)
        for i, ix in enumerate(sorted_indices):
            self.mses[i], self.eqs[i] = new_mse[ix], new_eqs[ix]

    def change(self, p_mute=0.05, p_rep=0.3):
        new_eqs, new_mse = copy.deepcopy(self.eqs), copy.deepcopy(self.mses)
        sorted_indices = np.argsort(new_mse)
        for i, ix in enumerate(sorted_indices):
            self.mses[i], self.eqs[i] = new_mse[ix], new_eqs[ix]

        for i in range(self.num):
            # 保留最好的那部分eqs不change，只cross over.
            if i < 1:
                continue

            # print(self.eqs[i].visualize())
            self.eqs[i].mutate(p_mute)
            replace_or_not = np.random.choice([False, True], p=([1 - p_rep, p_rep]))
            if replace_or_not:
                self.eqs[i].replace()
                self.mses[i], _ = evaluate_mse(self.eqs[i])


if __name__ == '__main__':
    # np.random.seed(10)
    # a_tree = Tree(max_depth=4, p_var=0.5)
    # print(is_an_equation(a_tree.preorder.split()))

    # pdb.set_trace()

    # pde = PDE(depth=4, max_width=3, p_var=0.5, p_mute=0.1)
    # evaluate_mse(pde)

    # pdb.set_trace()

    sga = SGA(num=20, depth=4, width=4, p_var=0.5, p_rep=1, p_mute=0.2, p_cro=0.5)
    sga.run(100)
