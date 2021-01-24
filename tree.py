from setup import *
import copy


class Node:
    """
        1. depth: 节点深度
        2. idx: 当前深度的第几个节点
        3. parent_idx: 父节点是上一层的第几个节点
        4. name: 节点详情
        5. child_num: 节点拥有几个孩子节点(单运算还是双运算)
        6. child_st: 下一层从第几个节点开始是当前的孩子节点
        7. var: 节点variable/operation
        8. cache: 保留运算至今的数据
        9. status: 初始化为child_num，用于记录遍历状态
        10. full: 完整信息，以OP或VAR形式表示
    """
    def __init__(self, depth, idx, parent_idx, name, full, child_num, child_st, var):
        self.depth = depth
        self.idx = idx
        self.parent_idx = parent_idx

        self.name = name
        self.child_num = child_num
        self.child_st = child_st
        self.status = self.child_num
        self.var = var
        self.full = full
        self.cache = copy.deepcopy(var)

    def __str__(self):
        return self.name

    def reset_status(self):
        self.status = self.child_num


class Tree: #对应于pde中的一个term
    def __init__(self, max_depth, p_var):
        self.max_depth = max_depth
        self.tree = [[] for i in range(max_depth)]
        self.preorder, self.inorder = None, None

        root = OPS[np.random.randint(0, len(OPS))] # e.g. ['sin', 1, np.sin], ['*', 2, np.multiply]
        node = Node(depth=0, idx=0, parent_idx=None, name=root[0], var=root[2], full=root,
                    child_num=int(root[1]), child_st=0)
        self.tree[0].append(node) # 初始节点

        depth = 1
        while depth < max_depth:
            next_cnt = 0 #child_st=next_cnt， child_st: 下一层从第几个节点开始是当前的孩子节点
            # 对应每一个父节点都要继续生成他们的子节点
            for parent_idx in range(len(self.tree[depth - 1])): #一个tree中某个depth处的node的子节点可以是多个，因此有可能在某个深度处存在多个node
                parent = self.tree[depth - 1][parent_idx] # 提取出对应深度处的对应操作符（某个node）
                if parent.child_num == 0: # 如果当前node没有子节点，则跳过当前循环的剩余语句，然后继续进行下一轮循环
                    continue
                for j in range(parent.child_num):
                    # parent var为d 且j为1时，必须确保右子节点为x
                    if parent.name == 'd' and j == 1: # j == 0 为d的左侧节点，j == 1为d的右侧节点
                        node = den[np.random.randint(0, len(den))] # 随机产生一个微分运算的denominator，一般是xyt
                        node = Node(depth=depth, idx=len(self.tree[depth]), parent_idx=parent_idx, name=node[0],
                                    var=node[2], full=node, child_num=int(node[1]), child_st=None)
                        self.tree[depth].append(node)
                    # 最底一层必须是var，不能是op
                    elif depth == max_depth - 1:
                        node = VARS[np.random.randint(0, len(VARS))]
                        node = Node(depth=depth, idx=len(self.tree[depth]), parent_idx=parent_idx, name=node[0],
                                    var=node[2], full=node, child_num=int(node[1]), child_st=None)
                        self.tree[depth].append(node)
                    else:
                        # 不是最底一层，p_var概率var， 区别在于child_st
                        if np.random.random() <= p_var:
                            node = VARS[np.random.randint(0, len(VARS))]
                            node = Node(depth=depth, idx=len(self.tree[depth]), parent_idx=parent_idx, name=node[0],
                                        var=node[2], full=node, child_num=int(node[1]), child_st=None)
                            self.tree[depth].append(node)
                        else:
                            node = OPS[np.random.randint(0, len(OPS))]
                            node = Node(depth=depth, idx=len(self.tree[depth]), parent_idx=parent_idx, name=node[0],
                                        var=node[2], full=node, child_num=int(node[1]), child_st=next_cnt)
                            next_cnt += node.child_num
                            self.tree[depth].append(node)
            depth += 1

        ret = []
        dfs(ret, self.tree, depth=0, idx=0)
        self.preorder = ' '.join([x for x in ret])
        model_tree = copy.deepcopy(self.tree)
        self.inorder = tree2str_merge(model_tree)

        # print(self.preorder)
        # print(self.inorder)
        # print('---------------')

    def mutate(self, p_mute): #直接替换原有tree中的某个节点，用同类型节点替换，因此后续位置不需要重新生成（类似替换了一个基因，而不是把后续基因序列重新产生，具有物理含义，也易于实现）
        global see_tree
        see_tree = copy.deepcopy(self.tree)
        depth = 1
        while depth < self.max_depth:
            next_cnt = 0
            idx_this_depth = 0  # 这个深度第几个节点
            for parent_idx in range(len(self.tree[depth - 1])):
                parent = self.tree[depth - 1][parent_idx]
                if parent.child_num == 0:
                    continue
                for j in range(parent.child_num):  # parent 的第j个子节点
                    # 不突变则跳过
                    not_mute = np.random.choice([True, False], p=([1 - p_mute, p_mute]))
                    if not_mute:
                        next_cnt += self.tree[depth][parent.child_st + j].child_num
                        continue

                    # 当前节点的类型
                    current = self.tree[depth][parent.child_st + j]
                    temp = self.tree[depth][parent.child_st + j].name
                    num_child = self.tree[depth][parent.child_st + j].child_num  # 当前变异节点的子节点数

                    # print('mutate!')

                    # 叶节点必须是var，不能是op
                    if num_child == 0:
                        node = VARS[np.random.randint(0, len(VARS))]
                        while node[0] == temp or (parent.name == 'd' and node[0] not in den[:, 0]):
                            if cheat and parent.name == 'd' and node[0] == 'x':
                                break
                            # 避免重复 或者 自己是d的右子节点而自己又不是x
                            node = VARS[np.random.randint(0, len(VARS))] # 重新抽取一个vars
                        new_node = Node(depth=depth, idx=idx_this_depth, parent_idx=parent_idx, name=node[0],
                                    var=node[2], full=node, child_num=int(node[1]), child_st=None)
                        self.tree[depth][parent.child_st + j] = new_node #替换成变异的节点
                    else:
                        if num_child == 1:
                            node = OP1[np.random.randint(0, len(OP1))]
                            while node[0] == temp:  # 避免重复
                                node = OP1[np.random.randint(0, len(OP1))]
                        elif num_child == 2:
                            node = OP2[np.random.randint(0, len(OP2))]
                            right = self.tree[depth + 1][current.child_st + 1].name
                            # 避免重复，避免生成d以打乱树结构（新d的右子节点不是x）
                            while node[0] == temp or (node[0] == 'd' and right not in den[:, 0]):
                                node = OP2[np.random.randint(0, len(OP2))]
                        else:
                            raise NotImplementedError("Error occurs!")

                        new_node = Node(depth=depth, idx=idx_this_depth, parent_idx=parent_idx, name=node[0],
                                    var=node[2], full=node, child_num=int(node[1]), child_st=next_cnt)
                        next_cnt += new_node.child_num
                        self.tree[depth][parent.child_st + j] = new_node

                    idx_this_depth += 1

            depth += 1

        ret = []
        dfs(ret, self.tree, depth=0, idx=0)
        self.preorder = ' '.join([x for x in ret])
        model_tree = copy.deepcopy(self.tree)
        self.inorder = tree2str_merge(model_tree)

        # print(self.preorder)
        # print(self.inorder)


def dfs(ret, a_tree, depth, idx): #辅助中序遍历，产生一个描述这个tree的名称序列
    # print(depth, idx)  # 深度优先遍历的顺序
    node = a_tree[depth][idx]
    ret.append(node.name) #记录当前操作
    for ix in range(node.child_num):
        if node.child_st is None:
            continue
        dfs(ret, a_tree, depth+1, node.child_st + ix) #进入下一层中下一个节点对应的子节点


def tree2str_merge(a_tree):
    for i in range(len(a_tree) - 1, 0, -1):
        for node in a_tree[i]:
            if node.status == 0:
                if a_tree[node.depth-1][node.parent_idx].status == 1:
                    a_tree[node.depth-1][node.parent_idx].name = a_tree[node.depth-1][node.parent_idx].name + ' ' + node.name + ')' \
                        if a_tree[node.depth-1][node.parent_idx].child_num == 2 \
                        else '( ' + a_tree[node.depth-1][node.parent_idx].name + ' ' + node.name + ')'
                elif a_tree[node.depth-1][node.parent_idx].status > 1:
                    a_tree[node.depth-1][node.parent_idx].name = '(' + node.name + ' ' + a_tree[node.depth-1][node.parent_idx].name
                else:
                    pdb.set_trace()
                    raise NotImplemented()
                a_tree[node.depth-1][node.parent_idx].status -= 1
    return a_tree[0][0].name


class Point:
    def __init__(self, idx, name, child_num, child_idx=[]):
        """
            1. idx: 当前序列的第几个节点
            2. parent_idx: 父节点是第几个节点
            3. name: 节点名称
            4. child_num: 节点拥有几个孩子节点
            5. child_idx: 孩子节点是序列的第几个
        """
        self.idx = idx
        self.name = name
        self.child_num = child_num
        self.child_idx = child_idx

    def __str__(self):
        return self.name

    def add_child(self, ix):
        self.child_idx.append(ix)


def is_an_equation(seq):  # e.g. (+ u - u u)
    def split(seq, idx):
        # last element is an op
        if idx >= len(seq): return np.inf

        # idx is the current node
        op = ALL[:, 0]
        root = ALL[np.where(op == seq[idx])][0]
        node = Point(idx=idx, name=root[0], child_num=int(root[1]))

        if node.child_num != 0:
            node.child_idx.append(idx + 1)  # might be wrong for the last node, not fatal though
            new_idx = split(seq, idx + 1)  # first child
            if node.child_num != 1:  # other children
                node.child_idx.append(new_idx)
                new_idx = split(seq, new_idx)
            return new_idx

        return idx + 1

    idx = 0
    end_idx = split(seq, idx)
    if end_idx != len(seq):
        return False
    return True


if __name__ == '__main__':
    tree = Tree(max_depth=4, p_var=0.5)
    print(tree.inorder)
    tree.mutate(p_mute=1)
    print(tree.inorder)

