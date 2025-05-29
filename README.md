# CLRS-lab2
## 算法伪代码
~~~ python
function main():
    read query_input
    read reference_input
    create GenomeMatcher with query_input and reference_input
    call align()

class GenomeMatcher:
    constants:
        MAX_JUMP = 1000
        MIN_SEG_LEN = 15

    constructor(query, reference):
        store query as list
        store reference as list
        create reverse complement of reference
        initialize empty memoization table (node → cost and parent)
        initialize empty set of processed nodes
        initialize processing queue

    function align():
        initialize starting CONTROL node at (0, 0)
        store it in memo with cost 0 and itself as parent
        enqueue the start node

        while the queue is not empty:
            pop current node from queue
            if query index of current node reaches end of query:
                reinsert node at front and break

            if node is already processed, skip it
            mark node as processed
            get cost from memo

            if node is CONTROL:
                for each reference index within MAX_JUMP range of query index:
                    enqueue FORWARD and REVERSE nodes starting from (ref_idx, qry_idx) with cost +1

            else if node is FORWARD:
                if both indices are valid and bases match:
                    enqueue next FORWARD node with no penalty and urgent priority
                else:
                    enqueue with penalty +1
                also try skipping one base in query or reference
                finally, allow returning to CONTROL mode with cost +1

            else if node is REVERSE:
                if reverse complement base matches query base:
                    enqueue next REVERSE node with no penalty and urgent priority
                else:
                    enqueue with penalty +1
                also try skipping one base in query or reference
                finally, allow returning to CONTROL mode with cost +1

        if queue is empty:
            print "[]"
            return

        set node to the front of the queue
        initialize empty segment list
        record current query and reference end positions

        while node is not the starting CONTROL(0, 0):
            get previous node from memo
            if switching from CONTROL to FORWARD or REVERSE:
                update r_end and q_end to current position
            if switching from FORWARD/REVERSE to CONTROL and segment is long enough:
                record segment (q_start, q_end, r_start, r_end)
            move to parent node

        print recorded segments in reverse order

class SegmentNode:
    attributes: kind (CONTROL/FORWARD/REVERSE), ref_idx, qry_idx
    define hash and equality for dictionary storage

class PathInfo:
    attributes: cost, parent
~~~
## 时空复杂度
Q = len(query)
R = len(reference)
时间复杂度：O(R × Q)
空间复杂度：O(R × Q)
