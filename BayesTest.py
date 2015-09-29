class Main(object):
    def __init__(self):
        return

class BayesCalculator(object):
    def __init__(self):
        return

class Node(object):
    def get_parent(self, parent_name): 
        for p in self._parents:
            if p.parent.name == parent_name:
                return p

    def get_parent_value(self, parent_name):
        for p in self._parents:
            if p.parent.name == parent_name:
                return p.parent_value

    def get_child(self, child_name):
        for c in self._children:
            if c.child.name == child_name:
                return c

    def get_child_value(self, child_name):
        for c in self._children:
            if c.child.name == child_name:
                return c.child_value
    
    def add_cond_prob(self, vars, prob):
        var_parents = [] 
        for p in self.parents:
            if p.parent.name in vars:
                var_parents.append([p, True])
            else:
                var_parents.append([p, False])

        self.cond_prob_table.append([var_parents, prob])

    def print_cond_prob_table(self):
        print "Conditional probability table for: ", self.name
        print "Probability of ", self.name , "being true given: "
        for entry in range(0, len(self.cond_prob_table)):
            for parents_prob in range(0, len(self.cond_prob_table[entry]) - 1):
                for var_val in range(0, len(self.cond_prob_table[entry][parents_prob])):
                    print self.cond_prob_table[entry][parents_prob][var_val][0].parent.name, "=", self.cond_prob_table[entry][parents_prob][var_val][1] ,
            print self.cond_prob_table[entry][len(self.cond_prob_table[entry]) - 1]
        return

    # only print the entry matching the specified values (if any)
    def print_cpt_entry(self, vars_vals):
        print 'Conditional probability for:'
        for var_val in vars_vals:
            print var_val , '=', vars_vals.get(var_val) ,
        print ':', self.get_cpt_entry_prob(vars_vals)
        return

    #return the whole entry in the Conditional Probability Table of the node matching the variables and specified values
    def get_cpt_entry(self, vars_vals):
        notmatch = False
        for entry in range(0, len(self.cond_prob_table)):
            for parents in range(0, len(self.cond_prob_table[entry]) - 1):
                for var_arc in range(0, len(self.cond_prob_table[entry][parents])):
                    if not ((vars_vals.has_key(self.cond_prob_table[entry][parents][var_arc][0].parent.name)) 
                            and (self.cond_prob_table[entry][parents][var_arc][1] == vars_vals.get(self.cond_prob_table[entry][parents][var_arc][0].parent.name))
                            ):
                        notmatch = True
                        break
                if notmatch == True:
                    notmatch = False
                    break
                else: 
                    return self.cond_prob_table[entry]#[len(self.cond_prob_table[entry]) - 1]
                return

    # just return the probability of the variables and specified values
    def get_cpt_entry_prob(self, vars_vals):
        cpt_entry = self.get_cpt_entry(vars_vals)
        return cpt_entry[len(cpt_entry) - 1]

    def p(self, event, value = True):
        """Return the conditional probability
        P(X=value | parents=parent_values), where parent_values
        are the values of parents in event. (event must assign each
        parent a value.)
        """
        if len(self.cond_prob_table) > 0:
            ptrue = self.get_cpt_entry_prob(event)
        else:
            ptrue = self.prior_prob

        if value == False:
            return 1 - ptrue
        return ptrue

    def __init__(self, name):
        self.name = name
        self.parents = []
        self.children = []
        self.prior_prob = 0
        self.cond_prob_table = []

class Arc(object):

    def __init__(self):

        self.parent = Node(None)
        self.child = Node(None)
        self.child_value = False
        self.parent_value = False

class BNet(object):
    def __init__(self):
        self.nodes = []
        self.arcs = []
        return

    def variable_values(self, variable_name):
        "Return the domain of variable_name"
        return [True, False]

    def variable_node(self, variable_name):
        for n in self.nodes:
            if n.name == variable_name:
                return n
        raise Exception("No variable with the name: %s exists" % variable_name)

    def add_arc(self, parent_node, child_node):
        newarc = Arc()
        for p in self.nodes:
            if p.name == parent_node:
                newarc.parent = p
                p.children.append(newarc)
            if p.name == child_node:
                newarc.child = p
                p.parents.append(newarc)
        self.arcs.append(newarc)
        
    def add_node(self, node_name, prior_prob=None):
        newnode = Node(node_name)
        newnode.prior_prob = prior_prob
        self.nodes.append(newnode)


class Enumerate(object):
    """Return the conditional probability distribution of query_variable
    given evidence evidence, from Bayesian Network bayes_network
    >>> Enumerate.ask('variable name', dict({otherVariableName:T, otherOtherVariableName:F}), bayesianNetwork).show_approx()
    """
    @staticmethod
    def ask(query_variable, evidence, bayes_network):
        assert query_variable not in evidence, "Query variable must be distinct from evidence"
        Q = ProbDist(query_variable)
        for xi in bayes_network.variable_values(query_variable):
            current_query_val = {query_variable: xi}
            Q[xi] = Enumerate.all(bayes_network.nodes, dict(evidence.items() + current_query_val.items()), bayes_network)
        return Q.normalize()

    """Return the sum of those entries in P(vars | e{others})
    consistent with e, where P is the joint distribution represented
    by bn, and e{others} means e restricted to bn's other variables
    (the ones other than vars). Parents must precede children in vars."""
    @staticmethod
    def all(variables, evidence, bayes_network):
        if not variables:
            return 1.0
        Y, rest = variables[0], variables[1:]
        Ynode = bayes_network.variable_node(Y.name)
        if Y.name in evidence:
            return Ynode.p(evidence, evidence[Y.name]) * Enumerate.all(rest, evidence, bayes_network)
        else:
            return sum(Ynode.p(evidence, y) * Enumerate.all(rest, dict(evidence.items() + dict({Y.name: y}).items()), bayes_network)
                       for y in bayes_network.variable_values(Y.name))

class ProbDist(object):
    """A discrete probability distribution.  You name the random variable
    in the constructor, then assign and query probability of values.
    >>> P = ProbDist('Flip'); P['H'], P['T'] = 0.25, 0.75; P['H']
    0.25
    >>> P = ProbDist('X', {'lo': 125, 'med': 375, 'hi': 500})
    >>> P['lo'], P['med'], P['hi']
    (0.125, 0.375, 0.5)
    """
    def __init__(self, varname='?', freqs=None):
        """If freqs is given, it is a dictionary of value: frequency pairs,
        and the ProbDist then is normalized."""
        self.prob={} 
        self.varname=varname
        self.values=[]
        if freqs:
            for (v, p) in freqs.items():
                self[v] = p
            self.normalize()

    def __getitem__(self, val):
        "Given a value, return P(value)."
        try: return self.prob[val]
        except KeyError: return 0

    def __setitem__(self, val, p):
        "Set P(val) = p."
        if val not in self.values:
            self.values.append(val)
        self.prob[val] = p

    def normalize(self):
        """Make sure the probabilities of all values sum to 1.
        Returns the normalized distribution.
        Raises a ZeroDivisionError if the sum of the values is 0.
        >>> P = ProbDist('Flip'); P['H'], P['T'] = 35, 65
        >>> P = P.normalize()
        >>> print '%5.3f %5.3f' % (P.prob['H'], P.prob['T'])
        0.350 0.650
        """
        total = float(sum(self.prob.values()))
        if not (1.0-epsilon < total < 1.0+epsilon):
            for val in self.prob:
                self.prob[val] /= total
        return self

    def show_approx(self, numfmt='%.3g'):
        """Show the probabilities rounded and sorted by key, for the
        sake of portable doctests."""
        return ', '.join([('%s: ' + numfmt) % (v, p)
                          for (v, p) in sorted(self.prob.items())])

epsilon = 0.001

T, F = True, False

burglary = BNet()

burglary.add_node('Burglary', 0.001)
burglary.add_node('Earthquake', 0.002)
burglary.add_node('Alarm', 0)
burglary.add_node('JohnCalls', 0)
burglary.add_node('MaryCalls', 0)

burglary.add_arc('Burglary', 'Alarm')
burglary.add_arc('Earthquake', 'Alarm')
burglary.variable_node('Alarm').add_cond_prob(['Burglary', 'Earthquake'], 0.95)
burglary.variable_node('Alarm').add_cond_prob(['Burglary'], 0.94)
burglary.variable_node('Alarm').add_cond_prob(['Earthquake'], 0.29)
burglary.variable_node('Alarm').add_cond_prob([], 0.001) 

burglary.add_arc('Alarm', 'JohnCalls')
burglary.variable_node('JohnCalls').add_cond_prob(['Alarm'], 0.90)
burglary.variable_node('JohnCalls').add_cond_prob([], 0.05)

burglary.add_arc('Alarm', 'MaryCalls')
burglary.variable_node('MaryCalls').add_cond_prob(['Alarm'], 0.70)
burglary.variable_node('MaryCalls').add_cond_prob([], 0.01)

for node in burglary.nodes:
    node.print_cond_prob_table()

burglary.variable_node("Alarm").print_cpt_entry({'Burglary' : T, 'Earthquake' : F})

print 'Probability of Burglary given JohnCalls = T and MaryCalls = T'
print Enumerate.ask('Alarm', dict({'JohnCalls': T, 'MaryCalls':F}), burglary).show_approx()





        
