/*
 * Copyright © 2009 Ian D. Romanick
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef __glu3_h__
#define __glu3_h__

#include <GL/gl.h>
#include <GL/glext.h>

#define GLU3_VERSION_0_1

struct GLUmat4;

struct GLUvec4 {
        GLfloat values[4];

#ifdef __cplusplus

        inline GLUvec4(void)
        {
        }

        inline GLUvec4(GLfloat v)
        {
                values[0] = v;
                values[1] = v;
                values[2] = v;
                values[3] = v;
        }

        inline GLUvec4(GLfloat x , GLfloat y, GLfloat z, GLfloat w)
        {
                values[0] = x;
                values[1] = y;
                values[2] = z;
                values[3] = w;
        }

        inline GLUvec4(const GLUvec4 &v)
        {
                values[0] = v.values[0];
                values[1] = v.values[1];
                values[2] = v.values[2];
                values[3] = v.values[3];
        }

        GLUvec4 operator *(const GLUmat4 &) const;

        GLUvec4 operator *(const GLUvec4 &) const;

        GLUvec4 operator *(GLfloat) const;

        GLUvec4 operator +(const GLUvec4 &) const;

        GLUvec4 operator -(const GLUvec4 &) const;
#endif /* __cplusplus */
};


#ifdef __cplusplus
inline GLUvec4 operator *(GLfloat f, const GLUvec4 &v)
{
        return v * f;
}

inline GLUvec4 &operator +=(GLUvec4 &l, const GLUvec4 &r)
{
        l = l + r;
        return l;
}

inline GLUvec4 &operator -=(GLUvec4 &l, const GLUvec4 &r)
{
        l = l - r;
        return l;
}

inline GLUvec4 &operator *=(GLUvec4 &l, const GLUvec4 &r)
{
        l = l * r;
        return l;
}

inline GLUvec4 &operator *=(GLUvec4 &l, GLfloat r)
{
        l = l * r;
        return l;
}
#endif /* __cplusplus */


struct GLUmat4 {
        struct GLUvec4 col[4];

#ifdef __cplusplus

        inline GLUmat4(void)
        {
        }

        inline GLUmat4(const GLUvec4 & c0, const GLUvec4 & c1,
                       const GLUvec4 & c2, const GLUvec4 & c3)
        {
                col[0] = c0;
                col[1] = c1;
                col[2] = c2;
                col[3] = c3;
        }

        inline GLUmat4(const GLUmat4 &m)
        {
                col[0] = m.col[0];
                col[1] = m.col[1];
                col[2] = m.col[2];
                col[3] = m.col[3];
        }


        GLUvec4 operator *(const GLUvec4 &) const;

        GLUmat4 operator *(const GLUmat4 &) const;

        GLUmat4 operator *(GLfloat) const;

        GLUmat4 operator +(const GLUmat4 &) const;

        GLUmat4 operator -(const GLUmat4 &) const;
#endif  /* __cplusplus */
};

#define GLU_MAX_STACK_DEPTH 32

struct GLUmat4Stack {
        struct GLUmat4 stack[GLU_MAX_STACK_DEPTH];
        unsigned top;

#ifdef __cplusplus
        GLUmat4Stack() : top(0)
        {
                /* empty */
        }
#endif  /* __cplusplus */
};


struct GLUarcball {
        unsigned viewport_x;
        unsigned viewport_y;
        unsigned viewport_width;
        unsigned viewport_height;
        unsigned click_x;
        unsigned click_y;
#ifdef __cplusplus
        void viewport(unsigned x, unsigned y, unsigned width, unsigned height)
        {
                viewport_x = x;
                viewport_y = y;
                viewport_width = width;
                viewport_height = height;
        }

        void click(unsigned x, unsigned y)
        {
                click_x = x;
                click_y = y;
        }

        GLUmat4 drag(unsigned end_x, unsigned end_y);
#endif  /* __cplusplus */
};


#ifdef __cplusplus

class GLUshapeConsumer {
public:
        virtual void vertex(const GLUvec4 &position,
                            const GLUvec4 &normal,
                            const GLUvec4 &tangent,
                            const GLUvec4 &uv) = 0;

        virtual void begin_primitive(GLenum mode) = 0;

        virtual void index(unsigned idx) = 0;

        virtual void end_primitive(void) = 0;
};


class GLUshapeProducer {
public:
        virtual ~GLUshapeProducer()
        {
        }

        void orientation(bool outside);

        virtual unsigned vertex_count(void) const = 0;

        virtual unsigned element_count(void) const = 0;

        virtual unsigned primitive_count(void) const = 0;

        virtual void generate(GLUshapeConsumer *consumer) const = 0;

protected:
        GLUshapeProducer(void) :
          normals_point_out(true)
        {
        }

        bool normals_point_out;
};


class GLUsphereProducer : public GLUshapeProducer {
public:
        GLUsphereProducer(GLdouble radius, GLint slices, GLint stacks);
        virtual unsigned vertex_count(void) const;
        virtual unsigned element_count(void) const;
        virtual unsigned primitive_count(void) const;
        virtual void generate(GLUshapeConsumer *consumer) const;

private:
        double radius;
        unsigned slices;
        unsigned stacks;
};


class GLUcubeProducer : public GLUshapeProducer {
public:
        GLUcubeProducer(GLdouble radius);
        virtual unsigned vertex_count(void) const;
        virtual unsigned element_count(void) const;
        virtual unsigned primitive_count(void) const;
        virtual void generate(GLUshapeConsumer *consumer) const;

private:
        double radius;
};
#endif  

#ifndef __cplusplus
typedef struct GLUvec4 GLUvec4;
typedef struct GLUmat4 GLUmat4;
typedef struct GLUmat4Stack GLUmat4Stack;
typedef struct GLUarcball GLUarcball;
#endif /*  __cplusplus */


#if defined(__cplusplus)
extern "C" {
#endif

GLfloat gluDot4_4v(const GLUvec4 *, const GLUvec4 *);

GLfloat gluDot3_4v(const GLUvec4 *, const GLUvec4 *);

GLfloat gluDot2_4v(const GLUvec4 *, const GLUvec4 *);

void gluCross4v(GLUvec4 *result, const GLUvec4 *u, const GLUvec4 *v);

void gluNormalize4v(GLUvec4 *result, const GLUvec4 *u);

GLfloat gluLength4v(const GLUvec4 *u);

GLfloat gluLengthSqr4v(const GLUvec4 *);

void gluOuter4v(GLUmat4 *result, const GLUvec4 *u, const GLUvec4 *v);


void gluMult4v_4v(GLUvec4 *result, const GLUvec4 *, const GLUvec4 *);

void gluDiv4v_4v(GLUvec4 *result, const GLUvec4 *, const GLUvec4 *);

void gluAdd4v_4v(GLUvec4 *result, const GLUvec4 *, const GLUvec4 *);

void gluSub4v_4v(GLUvec4 *result, const GLUvec4 *, const GLUvec4 *);

void gluMult4v_f(GLUvec4 *result, const GLUvec4 *, GLfloat);

void gluDiv4v_f(GLUvec4 *result, const GLUvec4 *, GLfloat);

void gluAdd4v_f(GLUvec4 *result, const GLUvec4 *, GLfloat);

void gluSub4v_f(GLUvec4 *result, const GLUvec4 *, GLfloat);

void gluMult4m_4m(GLUmat4 *result, const GLUmat4 *, const GLUmat4 *);

void gluAdd4m_4m(GLUmat4 *result, const GLUmat4 *, const GLUmat4 *);

void gluSub4m_4m(GLUmat4 *result, const GLUmat4 *, const GLUmat4 *);

void gluMult4m_4v(GLUvec4 *result, const GLUmat4 *m, const GLUvec4 *v);

void gluMult4m_f(GLUmat4 *result, const GLUmat4 *, GLfloat);

void gluScale4v(GLUmat4 *result, const GLUvec4 *u);

void gluTranslate3f(GLUmat4 *result, GLfloat x, GLfloat y, GLfloat z);

void gluTranslate4v(GLUmat4 *result, const GLUvec4 *v);
void gluRotate4v(GLUmat4 *result, const GLUvec4 *axis, GLfloat angle);

void gluLookAt4v(GLUmat4 *result, const GLUvec4 *eye, const GLUvec4 *center,
                 const GLUvec4 *up);

void gluFrustum6f(GLUmat4 *result, GLfloat left, GLfloat right, GLfloat bottom,
                  GLfloat top, GLfloat near, GLfloat far);

void gluPerspective4f(GLUmat4 *result, GLfloat fovy, GLfloat aspect,
                      GLfloat near, GLfloat far);

void gluOrtho4f(GLUmat4 *result, GLfloat left, GLfloat right, GLfloat bottom,
                GLfloat top);

void gluOrtho6f(GLUmat4 *result, GLfloat left, GLfloat right, GLfloat bottom,
                GLfloat top, GLfloat near, GLfloat far);
void gluTranspose4m(GLUmat4 *result, const GLUmat4 *m);

GLfloat gluDeterminant4_4m(const GLUmat4 *m);

GLboolean gluInverse4_4m(GLUmat4 *result, const GLUmat4 *m);


extern const GLUmat4 gluIdentityMatrix;

extern const GLchar *gluLoadTextFile(const char *file_name);

extern void gluUnloadTextFile(const GLchar *text);

extern void gluArcballViewport(GLUarcball *ball, unsigned x, unsigned y,
    unsigned width, unsigned height);

extern void gluArcballClick(GLUarcball *ball, unsigned start_x,
    unsigned start_y);

extern void gluArcballDrag(GLUarcball *ball, GLUmat4 *transformation,
    unsigned end_x, unsigned end_y);

#ifdef __cplusplus
};
#endif

#ifdef __cplusplus

GLfloat gluDot4(const GLUvec4 &, const GLUvec4 &);

GLfloat gluDot3(const GLUvec4 &, const GLUvec4 &);

GLfloat gluDot2(const GLUvec4 &, const GLUvec4 &);

inline GLUvec4 gluCross(const GLUvec4 &u, const GLUvec4 &v)
{
        GLUvec4 t;

        gluCross4v(& t, & u, & v);
        return t;
}

inline GLUvec4 gluNormalize(const GLUvec4 &v)
{
        GLUvec4 t;

        gluNormalize4v(& t, & v);
        return t;
}

inline GLfloat gluLength(const GLUvec4 &u)
{
        return gluLength4v(& u);
}

inline GLfloat gluLengthSqr(const GLUvec4 &u)
{
        return gluLengthSqr4v(& u);
}

inline GLUmat4 gluScale(const GLUvec4 &u)
{
        GLUmat4 result;

        gluScale4v(& result, & u);
        return result;
}

inline GLUmat4 gluScale(GLfloat x, GLfloat y, GLfloat z)
{
        GLUvec4 u(x, y, z, 1.0);
        GLUmat4 result;

        gluScale4v(& result, & u);
        return result;
}

inline GLUmat4 gluTranslate(GLfloat x, GLfloat y, GLfloat z)
{
        GLUmat4 result;

        gluTranslate3f(& result, x, y, z);
        return result;
}

inline GLUmat4 gluTranslate(const GLUvec4 &v)
{
        GLUmat4 result;

        gluTranslate4v(& result, & v);
        return result;
}
inline GLUmat4 gluRotate(const GLUvec4 &axis, GLfloat angle)
{
        GLUmat4 result;

        gluRotate4v(& result, & axis, angle);
        return result;
}

inline GLUmat4 gluLookAt(const GLUvec4 &eye, const GLUvec4 &center,
                         const GLUvec4 &up)
{
        GLUmat4 result;

        gluLookAt4v(& result, & eye, & center, & up);
        return result;
}

inline GLfloat gluDeterminant4(const GLUmat4 &m)
{
        return gluDeterminant4_4m(& m);
}

inline GLboolean gluInverse4(GLUmat4 &result, const GLUmat4 &m)
{
        return gluInverse4_4m(& result, & m);
}

inline GLUmat4 gluInverse4(const GLUmat4 &m)
{
        GLUmat4 result;

        gluInverse4_4m(& result, & m);
        return result;
}


inline GLUmat4 GLUarcball::drag(unsigned end_x, unsigned end_y)
{
        GLUmat4 result;

        gluArcballDrag(this, & result, end_x, end_y);
        return result;
}
#endif /* __cplusplus */

//#include "glu3_scalar.h"

#endif /* __glu3_h__ */
