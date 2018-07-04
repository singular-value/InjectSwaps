import math
from enum import Enum


class MappingMethodology(Enum):
    NAIVE_ROW_MAJOR = 1  # e.g. ABCD; EFGH; IJKL; ...
    GRAPH_MAPPER = 2  # invokes graph_mapper library
    ALTERNATING_ROW_MAJOR = 3  # useful for Ising Model (linear chain mapped to 2D)


def get_instructions_with_swaps(qasmf_filename,
                                mapping_methodology=MappingMethodology.NAIVE_ROW_MAJOR,
                                output_filename=None):
    """Returns transformed qasmf based on the provided file, but only permitting neighboring CNOTs.

    To facilitate CNOTs between qubits that are not neighbors (i.e. directly North South East or
    West of each other), SWAP operations are added to move the control next to the target.

    The returned qasmf specifies the physical coordinates of the qubits as opposed to their
    names originating from the Scaffold.

    The mapping methodology specifies the methodology for the initial mapping for the qubits.

    If an output filename is specified, the transformed qasmf will be written to the file.

    """
    physical_instructions = ""
    qubit_to_point = _get_initial_qubit_to_point(qasmf_filename, mapping_methodology)

    with open(qasmf_filename, 'r') as fp:
        for line in fp:
            if line.startswith('cbit'):
                physical_instructions += line
            elif _is_unary(line):
                physical_instructions += _get_transformed_unary_instruction(line, qubit_to_point)
            else:  # CNOT
                physical_instructions += _get_transformed_CNOT_instruction(line, qubit_to_point)

    if output_filename is not None:
        with open(output_filename, 'w') as fp:
            fp.write(physical_instructions)

    return physical_instructions


def _get_initial_qubit_to_point(qasmf_filename, mapping_methodology):
    qubits = _read_qubits(qasmf_filename)
    dim = math.ceil(len(qubits) ** 0.5)
    qubit_to_point = {}

    if mapping_methodology == MappingMethodology.NAIVE_ROW_MAJOR:
        x, y = 0, 0
        for qubit in qubits:
            qubit_to_point[qubit] = Point(x, y)
            x = (x + 1) % dim
            y = y + 1 if x == 0 else y
    elif mapping_methodology == MappingMethodology.GRAPH_MAPPER:  # not implemented yet
        raise NotImplementedError()
    elif mapping_methodology == MappingMethodology.ALTERNATING_ROW_MAJOR:
        x, y = 0, 0
        for qubit in qubits:
            qubit_to_point[qubit] = Point(x, y)
            # for even rows, go right. For odd rows, go left.
            if y % 2 == 0:
                x = (x + 1) % dim
                if x == 0:
                    x, y = dim - 1, y + 1
            else:
                x = (x - 1) % dim
                if x == dim - 1:
                    x, y = 0, y + 1

    return qubit_to_point


def _read_qubits(qasmf_filename):
    qubits = []
    with open(qasmf_filename, 'r') as fp:
        for line in fp:
            if _is_unary(line):
                if line[:2] in ['Rx', 'Ry', 'Rz']:
                    qubit = line.split()[1].split(',')[0]
                else:
                    qubit = line.split()[1]
                if qubit not in qubits:
                    qubits.append(qubit)
            else:  # CNOT
                qubit1, qubit2 = line.split()[1].split(',')
                if qubit1 not in qubits:
                    qubits.append(qubit1)
                if qubit2 not in qubits:
                    qubits.append(qubit2)
    return qubits


def _is_unary(instruction):
    return not instruction.startswith('CNOT')


def _get_transformed_unary_instruction(line, qubit_to_point):
    # rename the qubit in instruction to qubit_x_y:
    if line[:2] in ['Rx', 'Ry', 'Rz']:
        operator, (qubit, angle) = line.split()[0], line.split()[1].split(',')
        return '%s %s,%s\n' % (operator, qubit_to_point[qubit].coordinate_qubit(), angle)
    else:
        operator, qubit = line.split()
        return '%s %s\n' % (operator, qubit_to_point[qubit].coordinate_qubit())


def _get_transformed_CNOT_instruction(line, qubit_to_point):
    physical_instructions = ""
    operator, (control, target) = line.split()[0], line.split()[1].split(',')

    while(not _is_manhattan_neighbor(control, target, qubit_to_point)):
        physical_instructions += _swap_control_towards_target(control, target, qubit_to_point)

    physical_instructions += '%s %s,%s\n' % (operator, qubit_to_point[control].coordinate_qubit(),
                                           qubit_to_point[target].coordinate_qubit())

    return physical_instructions


def _is_manhattan_neighbor(qubit1, qubit2, qubit_to_point):
    """Returns True iff is qubit1 is one unit North, South, East, or West of qubit2."""
    point1, point2 = qubit_to_point[qubit1], qubit_to_point[qubit2]
    return abs(point1.x - point2.x) + abs(point1.y - point2.y) == 1


def _swap_control_towards_target(control, target, qubit_to_point):
    """Returns SWAP instruction to move control one unit towards target and updates qubit_to_point.
    
    Precondition: control and target must not be manhattan neighbors.

    """
    control_point, target_point = qubit_to_point[control], qubit_to_point[target]

    # perform swap in the direction of greater imbalance:
    if abs(control_point.x - target_point.x) >= abs(control_point.y - target_point.y):
        delta = -1 if control_point.x > target_point.x else 1
        neighbor_point = Point(control_point.x + delta, control_point.y)
    else:
        delta = -1 if control_point.y > target_point.y else 1
        neighbor_point = Point(control_point.x, control_point.y + delta)

    point_to_qubit = {v: k for k, v in qubit_to_point.items()}  # inefficient but fine
    neighbor = point_to_qubit[neighbor_point]

    qubit_to_point[control] = neighbor_point 
    qubit_to_point[neighbor] = control_point
    return 'SWAP %s,%s\n' % (control_point.coordinate_qubit(), neighbor_point.coordinate_qubit())


class Point(object):
    """Represents a geometric point in 2-D space.

    Direction of coordinate system is as follows (0,0 is top left):
        0,0  1,0  2,0
        0,1  1,1  2,1
        0,2  1,2  2,2
    """
    def __init__(self, x, y):
        assert isinstance(x, int), 'x (%s) needs to be an int' % x
        assert isinstance(y, int), 'y (%s) needs to be an int' % y
        self.x = x
        self.y = y

    def __repr__(self):
        return '(%s, %s)' % (self.x, self.y)

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, another):
        return self.x == another.x and self.y == another.y

    def coordinate_qubit(self):
        return 'qubit_%s_%s' % (self.x, self.y)
