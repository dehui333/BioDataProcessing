#include <algorithm>
#include <assert.h>
#include <Python.h>
#include <iostream>
#include <string>
#include <vector>

#include "models.h"

#define TYPE_COV 0 // no/low coverage
#define TYPE_INS 1
#define TYPE_DEL 2
#define TYPE_MIS 4
#define TYPE_HP 5
#define TYPE_NONE 6
#define TYPE_CLEAR 7
#define A_INT 1
#define C_INT 2
#define G_INT 4
#define T_INT 8
#define N_INT 15

/*

- Try to have separate options for skipping HP ins/HP del
- Try to change to smaller datatypes for the parsetuple if have time or sth
- I don't think the assert statements check anything. Try with print out or sth.
*/
extern uint8_t min_mapping_quality_global;

inline bool is_HP(const char *contig, std::uint32_t pos, std::uint32_t contig_len)
{
    std::uint8_t i;

    bool hp = true;

    if (pos + 3 >= contig_len)
    {
        hp = false;
    }
    else
    {
        for (i = 2; i <= 3; i++)
        {
            if (contig[pos + 1] != contig[pos + i])
            {
                hp = false;
            }
        }
    }
    if (hp)
        return true;

    hp = true;
    if (pos < 3)
    {
        hp = false;
    }
    else
    {
        for (i = 2; i <= 3; i++)
        {
            if (contig[pos - 1] != contig[pos - i])
            {

                hp = false;
            }
        }
    }
    return hp;
}

static PyObject *get_diff_with_assm_cpp(PyObject *self, PyObject *args)
{
    const char *bam_path;
    const char *contig_name;
    const char *contig;
    long contig_len; // There's some problem with these and ParseTuple if I try to use unsigned..
    long min_mapq;   // ^ it's due to using wrong type for boolean. should just use int or sth
    long min_coverage;
    double min_diff_prop;
    int skip_HP;
    if (!PyArg_ParseTuple(args, "sssllldi", &bam_path, &contig_name, &contig,
                          &contig_len, &min_mapq, &min_coverage, &min_diff_prop, &skip_HP))
        return NULL;
    auto bam = readBAM(bam_path);
    std::string contig_name_string{contig_name};
    auto pileup_iter = bam->pileup(contig_name_string + ":");
    min_mapping_quality_global = min_mapq;
    struct diff
    {
        std::uint32_t pos;
        std::uint8_t type;
        std::uint32_t len; // does not indicate ins len but len of stretch of positions on contig
        explicit diff(std::uint32_t pos, std::uint8_t type, std::uint32_t len) : pos(pos), type(type), len(len) {}
        inline void extend(std::uint32_t len = 1) { this->len += len; }
    };
    std::vector<diff> diffs;
    auto extend = [&diffs](std::uint32_t pos, std::uint8_t type, std::uint32_t len)
    {
        if (!diffs.empty() && diffs.back().type == type && (diffs.back().pos + diffs.back().len) == pos)
        {
            diffs.back().extend(len);
        }
        else
        {
            diffs.emplace_back(pos, type, len);
        }
    };
    long expected_rpos = 0;

    // need to do sth about has_next but 0 diff
    if (!pileup_iter->has_next())
    {
        extend(0, TYPE_NONE, 0);
    }
    while (pileup_iter->has_next())
    {
        // std::cout << "BEFORE " << std::endl;
        auto column = pileup_iter->next();
        long rpos = column->position;
        if (rpos < pileup_iter->start())
            continue;
        if (rpos >= pileup_iter->end())
            break;
        // std::cout <<"PASS" << std::endl;

        if (rpos > expected_rpos)
        {
            extend(expected_rpos, TYPE_COV, rpos - expected_rpos);
        }
        expected_rpos = rpos + 1;
        std::uint16_t coverage = column->count();

        if (coverage < min_coverage)
        {
            extend(rpos, TYPE_COV, 1);
            continue;
        }
        // std::cout << rpos << ", " << coverage << std::endl;
        std::uint16_t num_del = 0;
        std::uint16_t num_ins = 0;
        std::uint16_t counts[16];
        counts[A_INT] = 0;
        counts[C_INT] = 0;
        counts[G_INT] = 0;
        counts[T_INT] = 0;
        counts[N_INT] = 0;
        char int2char[16] = {
            'X', 'A', 'C', 'X',
            'G', 'X', 'X', 'X',
            'T', 'X', 'X', 'X',
            'X', 'X', 'X', 'N'};
        while (column->has_next())
        {
            auto r = column->next();
            if (r->is_refskip())
                continue;
            if (r->is_del())
            {
                // DELETION
                num_del++;
            }
            else
            {
                // POSITION
                auto ibase = r->ibase(0);
                counts[ibase]++;

                if (r->indel() > 0)
                {
                    num_ins++;
                }
            }
        }

        assert(coverage == counts[A_INT] + counts[C_INT] + counts[G_INT] + counts[T_INT] + counts[N_INT] + num_del);
        if ((double)num_del / coverage >= min_diff_prop)
        {
            if (is_HP(contig, rpos, contig_len))
            {
                if (!skip_HP)
                    extend(rpos, TYPE_HP, 1);
            }
            else
            {
                extend(rpos, TYPE_INS, 1);
            }

            continue;
        }
        std::uint8_t idx[5] = {A_INT, C_INT, G_INT, T_INT, N_INT};
        bool has_mis = false;
        for (auto i = 0; i < 5; i++)
        {
            if (int2char[idx[i]] != contig[rpos] && (double)counts[idx[i]] / coverage >= min_diff_prop)
            {
                has_mis = true;
                break;
            }
        }
        if (has_mis)
        {
            extend(rpos, TYPE_MIS, 1);
            continue;
        }

        if ((double)num_ins / coverage >= min_diff_prop)
        {
            if (is_HP(contig, rpos, contig_len))
            {
                if (!skip_HP)
                    extend(rpos, TYPE_HP, 1);
            }
            else
            {
                extend(rpos, TYPE_DEL, 1);
            }

            continue;
        }
    }
    if (diffs.empty())
    {
        extend(0, TYPE_CLEAR, 0);
    }
    PyObject *result_list = PyList_New(diffs.size());
    for (std::uint32_t i = 0; i < diffs.size(); i++)
    {
        PyObject *diff_tuple = PyTuple_New(3);
        PyTuple_SetItem(diff_tuple, 0, PyLong_FromUnsignedLong(diffs[i].pos));
        PyTuple_SetItem(diff_tuple, 1, PyLong_FromUnsignedLong(diffs[i].type));
        PyTuple_SetItem(diff_tuple, 2, PyLong_FromUnsignedLong(diffs[i].len));
        PyList_SetItem(result_list, i, diff_tuple);
    }

    return result_list;
}
/*
I think shifting ins/del can lead to a position being included wrongly...
* but im skipping indel anyway due to deletions leading to high counts as
we are not counting the number of occurrences but number of positions.
*/
static PyObject *get_snp_pos_cpp(PyObject *self, PyObject *args)
{
    const char *bam_path;
    const char *contig_name;
    const char *contig;
    long contig_len;
    long min_mapq;
    long min_coverage;
    double min_alt_prop;
    int skip_HP;
    int skip_indel;
    // bool skip_HP;
    if (!PyArg_ParseTuple(args, "sssllldii", &bam_path, &contig_name, &contig,
                          &contig_len, &min_mapq, &min_coverage, &min_alt_prop, &skip_HP, &skip_indel))
        return NULL;
    auto bam = readBAM(bam_path);
    std::string contig_name_string{contig_name};
    auto pileup_iter = bam->pileup(contig_name_string + ":");
    min_mapping_quality_global = min_mapq;
    std::vector<std::uint32_t> positions;
    while (pileup_iter->has_next())
    {
        // std::cout << "BEFORE " << std::endl;
        auto column = pileup_iter->next();
        long rpos = column->position;
        if (rpos < pileup_iter->start())
            continue;
        if (rpos >= pileup_iter->end())
            break;
        // std::cout <<"PASS" << std::endl;

        std::uint16_t coverage = column->count();

        if (coverage < min_coverage)
        {
            continue;
        }
        // std::cout << rpos << ", " << coverage << std::endl;
        std::uint16_t num_del = 0;
        std::uint16_t num_ins = 0;
        std::uint16_t counts[16];
        counts[A_INT] = 0;
        counts[C_INT] = 0;
        counts[G_INT] = 0;
        counts[T_INT] = 0;
        counts[N_INT] = 0;
        /*char int2char[16] = {
            'X', 'A', 'C', 'X',
            'G', 'X', 'X', 'X',
            'T', 'X', 'X', 'X',
            'X', 'X', 'X', 'N'};*/
        while (column->has_next())
        {
            auto r = column->next();
            if (r->is_refskip())
                continue;
            if (r->is_del())
            {
                // DELETION
                //num_del++;
                // del counted in number of occurrences, not number of bases now
                // do nothing here to avoid repeat counting del segments.
            }
            else
            {
                // POSITION
                auto ibase = r->ibase(0);
                counts[ibase]++;

                if (!skip_indel)
                {
                    if (r->indel() > 0)
                    {
                        num_ins++;
                    }
                    else if (r->indel() < 0)
                    {
                        num_del++;
                    }
                    
                }
            }
        }
        assert(coverage == counts[A_INT] + counts[C_INT] + counts[G_INT] + counts[T_INT] + counts[N_INT] + num_del);

        std::vector<std::uint16_t> all_counts;
        all_counts.reserve(7);
        all_counts.push_back(counts[A_INT]);
        all_counts.push_back(counts[C_INT]);
        all_counts.push_back(counts[G_INT]);
        all_counts.push_back(counts[T_INT]);
        all_counts.push_back(counts[N_INT]);

        if (!skip_indel && (!skip_HP || !is_HP(contig, rpos, contig_len)))
        {
            all_counts.push_back(num_del);
            all_counts.push_back(num_ins);
        }
        struct sort_struct
        {
            bool operator()(std::uint16_t i, std::uint16_t j)
            {
                return i > j;
            }
        } descending;
        std::sort(all_counts.begin(), all_counts.end(), descending);
        std::uint32_t sum_counts = counts[A_INT] + counts[C_INT] + counts[G_INT] + counts[T_INT] + counts[N_INT] + num_del + num_ins;
        if ((double)all_counts[1] / sum_counts > min_alt_prop)
        {
            positions.push_back(rpos);
        }
    }

    PyObject *result_list = PyList_New(positions.size());
    for (std::uint32_t i = 0; i < positions.size(); i++)
    {
        PyList_SetItem(result_list, i, PyLong_FromUnsignedLong(positions[i]));
    }

    return result_list;

    // Py_RETURN_NONE;
}

static PyMethodDef check_alignment_methods[] = {
    {"get_diff_with_assm", get_diff_with_assm_cpp, METH_VARARGS, "Return apparent disagreements between reads and assembly."},
    {"get_snp_pos", get_snp_pos_cpp, METH_VARARGS, "Return apparent snps positions."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef check_alignment_definition = {
    PyModuleDef_HEAD_INIT,
    "check_alignment",
    "Check alignments.",
    -1,
    check_alignment_methods};

PyMODINIT_FUNC PyInit_check_alignment(void)
{
    Py_Initialize();
    return PyModule_Create(&check_alignment_definition);
}
