from dataclasses import dataclass
import gurobipy as gp
from gurobipy import GRB



@dataclass
class VertexProfile:  # vertex profile in the LP
    deg: int  # degree in EDCS
    is_matched: int  # 0-1 value, is it matched in the matching of EDCS
    is_opt: int # 0-1 value, is it matched in maximum matching

    # Hall's witness sets. L = A \cup Ab. R = NA \cup NAb. There is no edge of EDCS between A and NAb.
    # All vertices of NA and Ab should be matched
    part: str
    possible_parts = ['A', 'NA', 'Ab', 'NAb']


    # this function check if the constructed vertex profile is valid
    def is_valid(self) -> bool:
        if self.deg > beta or self.deg < 0:
            return False

        if self.is_matched not in [0, 1]:
            return False

        if self.part in ['NA', 'Ab'] and self.is_matched == 0:
            return False


        return True


@dataclass
class EdgeProfile:  # edge profile in the LP
    is_opt: int # 0-1 value, is it a maximum matching edge
    is_edcs: int # 0-1 value, is it an edge of EDCS
    is_matching: int # 0-1 value, is it the maximum matching of EDCS

    Lv: VertexProfile # left vertex profile of the edge, Lv = A \cup Ab
    Rv: VertexProfile # right vertex profile of the edge, Rv = NA \cup NAb


    # this function check if the constructed edge profile is valid
    def is_valid(self) -> bool:
        Lv, Rv = self.Lv, self.Rv
        if Lv.part not in ['A', 'Ab'] or Rv.part not in ['NA', 'NAb']:
            return False

        if (not self.is_opt) and (not self.is_edcs):
            return False

        if self.is_matching and (not self.is_edcs):
            return False

        if self.is_edcs and (Rv.deg + Lv.deg > beta):
            return False
        if (not self.is_edcs) and (Rv.deg + Lv.deg < beta_minus):
            return False

        if self.is_matching and ((not Rv.is_matched) or (not Lv.is_matched)):
            return False

        if self.is_edcs and self.Lv.part == 'A' and self.Rv.part == 'NAb':
            return False

        if self.is_matching and self.Rv.part == 'NA':
            if self.Lv.part != 'A':
                return False

        if self.is_matching and self.Lv.part == 'Ab':
            if self.Rv.part != 'NAb':
                return False

        if self.is_opt and ((not Rv.is_opt) or (not Lv.is_opt)):
            return False

        if self.is_edcs and ((Rv.deg == 0) or (Lv.deg == 0)):
            return False

        return True


def create_VertexProfiles():
    vertex_profiles = []
    for part in VertexProfile.possible_parts:
        for deg in range(beta):
            for is_matched in [0, 1]:
                for is_opt in [0, 1]:
                    vertex_profile = VertexProfile(deg, is_matched, is_opt, part)
                    if vertex_profile.is_valid():
                        vertex_profiles.append(vertex_profile)
    return vertex_profiles


def create_EdgeProfiles(vertex_profiles):
    global adj
    edge_profiles = []
    for Lv_ind, Lv in enumerate(vertex_profiles):
        for Rv_ind, Rv in enumerate(vertex_profiles):
            for is_opt in [0, 1]:
                for is_edcs in [0, 1]:
                    for is_matching in [0, 1]:
                        edge_profile = EdgeProfile(is_opt, is_edcs, is_matching, Lv, Rv)
                        if edge_profile.is_valid():
                            adj[Lv_ind].append((len(edge_profiles), edge_profile))
                            adj[Rv_ind].append((len(edge_profiles), edge_profile))
                            edge_profiles.append(edge_profile)

    return edge_profiles

for beta in range(2, 100):
    for beta_minus in range(1, beta):
        vertex_profiles = create_VertexProfiles()
        adj = [[] for i in range(len(vertex_profiles))]
        edge_profiles = create_EdgeProfiles(vertex_profiles)

        n1 = len(vertex_profiles)
        n2 = len(edge_profiles)
        n = n1 + n2

        constraints = []
        b = []

        # creating constraints of LP
        for vind, v in enumerate(vertex_profiles):
            matched = [0] * n
            edcs = [0] * n
            opt = [0] * n

            for eind, e in adj[vind]:
                if e.Lv == v or e.Rv == v:

                    if e.is_matching:
                        matched[n1 + eind] = 1

                    if e.is_edcs:
                        edcs[n1 + eind] = 1

                    if e.is_opt:
                        opt[n1 + eind] = 1

            if v.is_matched > 0:
                matched[vind] = -v.is_matched
                constraints.append(matched)

            if v.deg > 0:
                edcs[vind] = -v.deg
                constraints.append(edcs)

            if v.is_opt > 0:
                opt[vind] = -v.is_opt
                constraints.append(opt)


        b = [0] * len(constraints)

        constraint = [0] * n1
        for e in edge_profiles:
            constraint.append(e.is_matching)
        constraints.append(constraint)
        b.append(1)

        # setting the objective value of LP
        c = [0] * n1
        for e in edge_profiles:
            c.append(e.is_opt)

        num_variables = n
        num_constraints = len(b)

        model = gp.Model("EDCS LP")
        variables = model.addVars(range(n), lb=0.0)

        model.setObjective(gp.quicksum(c[i] * variables[i] for i in range(num_variables)), sense=GRB.MAXIMIZE)

        g_constraints = model.addConstrs((gp.quicksum(constraints[j][i] * variables[i] for i in range(num_variables)) == b[j] for j in range(num_constraints)))

        model.optimize()

        print(beta, beta_minus, 1./model.objVal)
