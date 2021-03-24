//
// Created by vfs on 2/5/2021.
//

// i couldn't figure out how to write to the SDL front buffer, i doubt it's easily doable since it's probably platform specific.
// so in order to have a constant aspect ratio i'll use a triple buffer.

// i give as input a vertex stream with positions in cartesian coordinates, normals and texture coordinates
// a valid vertex shader must return a position in homogenous coordinates, a normal and a texture coordinate
// a valid pixel shader must return a color


// includes

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <utility>

#include <SDL2/SDL.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "debugbreak.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

bool parse(const char *source, VertexList list);


// compiler definitions

//#define ASSERT(expr) if (!(expr)) *((int *) 0) = 0;
#define ARRAYCOUNT(arr) sizeof(arr) / sizeof(arr[0])

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;




// struct declarations

enum Primitive {
	POINT,
	LINE,
	TRIANGLE
};

struct Color {
	s16 r;
	s16 g;
	s16 b;
	s16 a;

	Color operator *(float k) const {
		return {(s16) ((float) r * k), (s16) ((float) g * k), (s16) ((float) b * k)};
	}

	Color operator -(Color c) const {
		return {(s16) (r - c.r), (s16) (g - c.g), (s16) (b - c.b)};
	}

	Color operator +(Color c) const {
		return {(s16) (r + c.r), (s16) (g + c.g), (s16) (b + c.b)};
	}
};

Color operator *(float k, Color c) {
	return {(s16) ((float) c.r * k), (s16) ((float) c.g * k), (s16) ((float) c.b * k)};
}

struct Vertex {
	glm::vec3 pos;
	glm::vec3 norm;
	glm::vec2 tex;
};

struct Vertex4 {
	glm::vec4 pos;
	glm::vec3 norm;
	glm::vec2 tex;
};

struct VertexList {
	Vertex *memory;
	u32 count;
	u32 size;
};

struct Vertex4List {
	Vertex4 *memory;
	u32 count;
	u32 size;
};

struct IndexList {
	u32 *memory;
	u32 count;
	u32 size;
};

struct Mesh {
	Vertex *vertices;
	u32 vertex_count;

	u32 *indices;
	u32 index_count;
};

struct ListMesh {
	VertexList vertexlist;
	IndexList indexlist;
};

struct Triangle {
	Vertex v0;
	Vertex v1;
	Vertex v2;
};

struct Transform {
	glm::vec3 translation;
	glm::vec3 rotation;
	glm::vec3 scale;
};

struct Camera {
	glm::vec3 translation;
	float pitch;
	float yaw;
};

struct Model {
	Mesh mesh;

	Transform t;
};

struct Plane {
	glm::vec3 n;
	float d;
};

struct Texture {
	u8 *data;
	s32 w;
	s32 h;
	s32 c;
};

struct Keyboard {
	bool forward;
	bool backward;
	bool left;
	bool right;
};


// forward fuction delcarations

void draw_indexed(Vertex *vertices, u32 vertex_count, u32 *indices, u32 index_count, Transform t);


// global memory

Color red = {255, 0, 0, 0};
Color green = {0, 255, 0, 0};
Color blue = {0, 0, 255, 0};
Color white = {255, 255, 255, 0};
Color black = {0, 0, 0, 0};
Color yellow = {255, 255, 0, 0};
Color purple = {255, 0, 255, 0};
Color cyan = {0, 255, 255, 0};


bool running = true;

SDL_bool rel_mouse = SDL_FALSE;

u64 frame;
u64 counter_freq;

float t = 0.0f;
const float dt = 1.0f / 240.0f;
float accumulator = 0.0f;

SDL_Window *window;
SDL_Surface *backbuffer;
SDL_Surface *triplebuffer;
float *depthbuffer;

Primitive mode;

s32 window_width = 1280;
s32 window_height = 720;

Camera camera;

float viewport_width = 2.0f;
float viewport_height = 2.0f;
float viewport_distance = 1.0f;

// should be even numbers
const s32 canvas_width = 500;
const s32 canvas_height = 500;

glm::mat4x3 projection_matrix;
glm::mat3x3 canvas_matrix;

Texture texture;

Keyboard keyboard;

Mesh cube_mesh;


// near plane has direction (0, 0, 1)
Plane near_plane = {glm::normalize(glm::vec3(0.0f, 0.0f, 1.0f)), -viewport_distance};
// far plane has direction (0, 0, -1)

// left plane has direction (1, 0, 1)
Plane left_plane = {glm::normalize(glm::vec3(1.0f, 0.0f, 1.0f)), 0.0f};
// right plane has direction (-1 , 0, 1)
Plane right_plane = {glm::normalize(glm::vec3(-1.0f, 0.0f, 1.0f)), 0.0f};

// top plane has direction (0, -1, 1)
Plane top_plane = {glm::normalize(glm::vec3(0.0f, -1.0f, 1.0f)), 0.0f};
// bottom plane has direction (0, 1, 1)
Plane bottom_plane = {glm::normalize(glm::vec3(0.0f, 1.0f, 1.0f)), 0.0f};




// function definitions

void camera_init(Camera *camera) {
	camera->translation = {0.0f, 0.0f, 0.0f};
	camera->pitch = glm::radians(0.0f);
	camera->yaw = glm::radians(0.0f);
}

glm::mat4x4 camera_get_matrix(Camera *camera) {


	glm::mat4x4 rotation;
	glm::mat4x4 translation;

	float rotation_array[16] = {
			cosf(camera->yaw), 0.0f, -sinf(camera->yaw), 0.0f,
			sinf(camera->yaw) * sinf(camera->pitch), cosf(camera->pitch), cosf(camera->yaw) * sinf(camera->pitch), 0.0f,
			sinf(camera->yaw) * cosf(camera->pitch), -sinf(camera->pitch), cosf(camera->yaw) * cos(camera->pitch), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
	};

	float translation_array[16] = {
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			camera->translation.x, camera->translation.y, camera->translation.z, 1.0f
	};

	rotation = glm::make_mat4x4(rotation_array);
	translation = glm::make_mat4x4(translation_array);

	rotation = glm::inverse(rotation);
	translation = glm::inverse(translation);

	glm::mat4x4 camera_matrix = rotation * translation;

	return camera_matrix;
}

void mesh_init(Mesh *mesh, Vertex *vertices, u32 vertex_count, u32 *indices, u32 index_count) {
	mesh->vertices = vertices;
	mesh->vertex_count = vertex_count;

	mesh->indices = indices;
	mesh->index_count = index_count;
}

void model_init(Model *model, Mesh mesh, Transform t) {
	model->mesh = mesh;

	model->t = t;
}

void model_render(Model *model) {
	draw_indexed(model->mesh.vertices, model->mesh.vertex_count, model->mesh.indices, model->mesh.index_count, model->t);
}

glm::mat4x4 transform_get_matrix(Transform *t) {
	glm::mat4x4 scale;
	glm::mat4x4 rotation;
	glm::mat4x4 translation;

	float scale_array[16] = {
			t->scale.x,       0.0f,       0.0f,       0.0f,
			0.0f,             t->scale.y, 0.0f,       0.0f,
			0.0f,             0.0f,       t->scale.z, 0.0f,
			0.0f,             0.0f,       0.0f,       1.0f
	};

	float rotation_array[16] = {
			cosf(t->rotation.z) * cosf(t->rotation.y), sinf(t->rotation.z) * cosf(t->rotation.y), -sinf(t->rotation.y), 0.0f,
			cosf(t->rotation.z) * sinf(t->rotation.y) * sinf(t->rotation.x) - sinf(t->rotation.z) * cosf(t->rotation.x), sinf(t->rotation.z) * sinf(t->rotation.y) * sinf(t->rotation.x) + cosf(t->rotation.z) * cosf(t->rotation.x), cosf(t->rotation.y) * sinf(t->rotation.x), 0.0f,
			cosf(t->rotation.z) * sinf(t->rotation.y) * cosf(t->rotation.x) + sinf(t->rotation.z) * sinf(t->rotation.x), sinf(t->rotation.z) * sinf(t->rotation.y) * cosf(t->rotation.x) - cosf(t->rotation.z) * sinf(t->rotation.x), cosf(t->rotation.y) * cosf(t->rotation.x), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
	};

	float translation_array[16] = {
			1.0f,             0.0f,             0.0f,             0.0f,
			0.0f,             1.0f,             0.0f,             0.0f,
			0.0f,             0.0f,             1.0f,             0.0f,
			t->translation.x, t->translation.y, t->translation.z, 1.0f
	};

	scale = glm::make_mat4x4(scale_array);
	rotation = glm::make_mat4x4(rotation_array);
	translation = glm::make_mat4x4(translation_array);

	glm::mat4x4 model_matrix = translation * rotation * scale;

	return model_matrix;
}

void vertexlist_init(VertexList *vertexlist) {
	vertexlist->memory = (Vertex *) malloc(8 * sizeof(Vertex));
	vertexlist->size = 8;
	vertexlist->count = 0;
}

s32 vertexlist_add(VertexList *vertexlist, Vertex vertex) {
	s32 ret = -1;

	if (vertexlist->count < vertexlist->size - 2) {
		ret = vertexlist->count;
		vertexlist->memory[vertexlist->count++] = vertex;
	} else {
		vertexlist->size *= 2;
		vertexlist->memory = (Vertex *) realloc(vertexlist->memory, vertexlist->size * sizeof(Vertex));
		ret = vertexlist->count;
		vertexlist->memory[vertexlist->count++] = vertex;
	}

	return ret;
}

Vertex vertexlist_get(VertexList *vertexlist, u32 index) {
	if (index < vertexlist->count) {
		return vertexlist->memory[index];
	} else {
		return {};
	}
}

void vertexlist_set(VertexList *vertexlist, u32 index, Vertex vertex) {
	if (index < vertexlist->count) {
		vertexlist->memory[index] = vertex;
	}
}

VertexList vertexlist_clone(VertexList *vertexlist) {
	VertexList ret = {};

	for (u32 i = 0; i < vertexlist->count; ++i) {
		vertexlist_add(&ret, vertexlist_get(vertexlist, i));
	}

	return ret;
}

void vertexlist_free(VertexList *vertexlist) {
	free(vertexlist->memory);
}

void vertex4list_init(Vertex4List *vertexlist) {
	vertexlist->memory = (Vertex4 *) malloc(8 * sizeof(Vertex4));
	vertexlist->size = 8;
	vertexlist->count = 0;
}

s32 vertex4list_add(Vertex4List *vertexlist, Vertex4 vertex) {
	s32 ret = -1;

	if (vertexlist->count < vertexlist->size - 2) {
		ret = vertexlist->count;
		vertexlist->memory[vertexlist->count++] = vertex;
	} else {
		vertexlist->size *= 2;
		vertexlist->memory = (Vertex4 *) realloc(vertexlist->memory, vertexlist->size * sizeof(Vertex4));
		ret = vertexlist->count;
		vertexlist->memory[vertexlist->count++] = vertex;
	}

	return ret;
}

Vertex4 vertex4list_get(Vertex4List *vertexlist, u32 index) {
	if (index < vertexlist->count) {
		return vertexlist->memory[index];
	} else {
		return {};
	}
}

void vertex4list_set(Vertex4List *vertexlist, u32 index, Vertex4 vertex) {
	if (index < vertexlist->count) {
		vertexlist->memory[index] = vertex;
	}
}

Vertex4List vertex4list_clone(Vertex4List *vertexlist) {
	Vertex4List ret = {};

	for (u32 i = 0; i < vertexlist->count; ++i) {
		vertex4list_add(&ret, vertex4list_get(vertexlist, i));
	}

	return ret;
}

void vertex4list_free(Vertex4List *vertexlist) {
	free(vertexlist->memory);
}

glm::vec3 homogenous_to_cartesian(glm::vec4 h) {
	if (h.w == 0.0f) return {h.x, h.y, h.z};

	return {h.x / h.w, h.y / h.w, h.z / h.w};
}

glm::vec4 cartesian_to_homogenous(glm::vec3 c) {
	return {c.x, c.y, c.z, 1.0f};
}

float signed_distance(Plane plane, glm::vec4 p) {
	glm::vec3 p0 = {p.x, p.y, p.z};

	float result = glm::dot(plane.n, p0) + plane.d;

	return result;
}

void update_backbuffer() {

	float triplebuffer_ar = (float) triplebuffer->w / (float) triplebuffer->h;
	float window_ar = (float) backbuffer->w / (float) backbuffer->h;

	if (triplebuffer_ar > window_ar) {
		s32 window_height_without_bars = (s32) ((float) backbuffer->w / (float) triplebuffer_ar);
		s32 bar_height = ((s32) window_height - window_height_without_bars) / 2;

		SDL_FillRect(backbuffer, NULL, 0x0);

		SDL_Rect destination = {0, bar_height, backbuffer->w, window_height_without_bars};
		SDL_BlitScaled(triplebuffer, NULL, backbuffer, &destination);
	} else if (triplebuffer_ar < window_ar) {
		s32 window_width_without_bars = (s32) (window_height * triplebuffer_ar);
		s32 bar_width = (s32) ((window_width - window_width_without_bars) / 2);

		SDL_FillRect(backbuffer, NULL, 0x0);

		SDL_Rect destination = {bar_width, 0, window_width_without_bars, backbuffer->h};
		SDL_BlitScaled(triplebuffer, NULL, backbuffer, &destination);
	} else {
		SDL_BlitScaled(triplebuffer, NULL, backbuffer, NULL);
	}
}

void set_pixel(int x, int y, float z, Color color) {
	assert(x >= 0);
	assert(x < triplebuffer->w);
	assert(y >= 0);
	assert(y < triplebuffer->h);

	y = (canvas_height - 1) - y;

	float *pixel_z = depthbuffer + (y * canvas_width + x);

	if (z > *pixel_z) {
		u32 color_as_u32 = (color.r << 16) | (color.g << 8) | color.b;

		u32 *memory = (u32 *) triplebuffer->pixels;
		u32 *pixel = memory + (y * triplebuffer->w + x);

		*pixel = color_as_u32;
		*pixel_z = z;
	}
}

void set_pixel_center(int x, int y, float z, Color color) {
	--x;
	--y;

	assert(x >= -triplebuffer->w / 2);
	assert(x < triplebuffer->w / 2);
	assert(y >= -triplebuffer->h / 2);
	assert(y < triplebuffer->h / 2);

	x += triplebuffer->w / 2;
	y += triplebuffer->h / 2;

	set_pixel(x, y, z, color);
}

void clear(Color color) {
	for (int x = 0; x < triplebuffer->w; ++x) {
		for (int y = 0; y < triplebuffer->h; ++y) {
			set_pixel(x, y, 0, color);
		}
	}
}

void clear_depthbuffer() {
	for (int x = 0; x < triplebuffer->w; ++x) {
		for (int y = 0; y < triplebuffer->h; ++y) {
			*(depthbuffer + (y * canvas_width + x)) = -INFINITY;
		}
	}
}

void draw_line(glm::ivec2 p0, glm::ivec2 p1, Color color) {
	bool steep = false;
	if (std::abs((int) (p1.y - p0.y)) > std::abs((int) (p1.x - p0.x))) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}

	if (p0.x > p1.x) {
		std::swap(p0.x, p1.x);
		std::swap(p0.y, p1.y);
	}

	int dx = p1.x - p0.x;
	int dy = p1.y - p0.y;

	int derror = std::abs(dy) * 2;
	int error = 0;

	int y = p0.y;

	if (steep) {
		for (int x = p0.x; x <= p1.x; ++x) {
			set_pixel_center(y, x, 1.0f, color);
			error += derror;
			if (error > dx) {
				y += (p1.y > p0.y ? 1 : -1);
				error -= dx * 2;
			}
		}
	} else {
		for (int x = p0.x; x <= p1.x; ++x) {
			set_pixel_center(x, y, 1.0f, color);
			error += derror;
			if (error > dx) {
				y += (p1.y > p0.y ? 1 : -1);
				error -= dx * 2;
			}
		}
	}
}

void draw_triangle_wireframe(Vertex p0, Vertex p1, Vertex p2) {
	draw_line(p0.pos, p1.pos, black);
	draw_line(p1.pos, p2.pos, white);
	draw_line(p2.pos, p0.pos, red);
}

Vertex4 clip_segment_against_plane(Plane plane, Vertex4 p0, Vertex4 p1) {
	float t = (-plane.d - glm::dot(plane.n, homogenous_to_cartesian(p0.pos))) / (glm::dot(plane.n, homogenous_to_cartesian(p1.pos) - homogenous_to_cartesian(p0.pos)));

	glm::vec4 q_pos = p0.pos + t * (p1.pos - p0.pos);
	glm::vec3 q_norm = p0.norm + t * (p1.norm - p0.norm);
	glm::vec2 q_tex = p0.tex + t * (p1.tex - p0.tex);

	return {q_pos, q_norm, q_tex};
}

Vertex4List clip_vertexlist_against_plane(Plane plane, Vertex4List vertexlist) {

	Vertex4List ret = {};
	vertex4list_init(&ret);

	if (vertexlist.count == 0) return ret;

	assert(vertexlist.count % 3 == 0);

	for (int i = 0; i < vertexlist.count - 2; i += 3) {
		Vertex4 v0 = vertex4list_get(&vertexlist, i);
		Vertex4 v1 = vertex4list_get(&vertexlist, i + 1);
		Vertex4 v2 = vertex4list_get(&vertexlist, i + 2);

		float d_0 = signed_distance(plane, v0.pos);
		float d_1 = signed_distance(plane, v1.pos);
		float d_2 = signed_distance(plane, v2.pos);


		if (d_0 < 0 && d_1 < 0 && d_2 < 0) {
			// nop
		} else if (d_0 >= 0 && d_1 >= 0 && d_2 >= 0) {
			vertex4list_add(&ret, v0);
			vertex4list_add(&ret, v1);
			vertex4list_add(&ret, v2);

		} else if (d_0 >= 0 && d_1 < 0 && d_2 < 0) {
			//intersect the segments p0->p1 and p0->p2

			// p0 -> p1
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v0, v1);

			// p0 -> p2
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v0, v2);

			vertex4list_add(&ret, v0);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, clipped_1);

		} else if (d_0 < 0 && d_1 >= 0 && d_2 < 0) {
			//intersect the segments p1->p0 and p1->p2

			// p1 -> p0
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v1, v0);


			// p1 -> p2
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v1, v2);

			vertex4list_add(&ret, v1);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, clipped_1);

		} else if (d_0 < 0 && d_1 < 0 && d_2 >= 0) {
			//intersect the segments p2->p0 and p2->p1

			// p2 -> p0
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v2, v0);

			// p2 -> p1
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v2, v1);

			vertex4list_add(&ret, v2);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, clipped_1);

		} else if (d_0 >= 0 && d_1 >= 0 && d_2 < 0) {
			//intersect the segments p0->p2 and p1->p2

			// p0 -> p2
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v0, v2);

			// p1 -> p2
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v1, v2);

			// p0 is A, p1 is B, p2 is C
			// p0 -> clipped_0 -> p1 is one triangle
			// p1 -> clipped_0 -> clipped_1 is the other

			vertex4list_add(&ret, v0);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, v1);
			vertex4list_add(&ret, v1);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, clipped_1);

		} else if (d_0 >= 0 && d_1 < 0 && d_2 >= 0) {
			//intersect the segments p0->p1 and p2->p1

			// p0 -> p1
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v0, v1);

			// p2 -> p1
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v2, v1);

			// p0 is B, p1 is C, p2 is A
			// p2 -> clipped_0 -> p0
			// p2 -> clipped_0 -> clipped_1

			vertex4list_add(&ret, v0);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, v2);
			vertex4list_add(&ret, v2);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, clipped_1);
		} else if (d_0 < 0 && d_1 >= 0 && d_2 >= 0) {
			//intersect the segments p1->p0 and p2->p0

			// p1 -> p0
			Vertex4 clipped_0 = clip_segment_against_plane(plane, v1, v0);

			// p2 -> p0
			Vertex4 clipped_1 = clip_segment_against_plane(plane, v2, v0);

			// p0 is C, p1 is A, p2 is B
			// p1 -> clipped_0 -> p2
			// p2 -> clipped_0 -> clipped_1

			vertex4list_add(&ret, v1);
			vertex4list_add(&ret, clipped_0);
			vertex4list_add(&ret, v2);
			vertex4list_add(&ret, v2);
			vertex4list_add(&ret, clipped_1);
			vertex4list_add(&ret, clipped_0);
		} else {
			printf("MAYBE A VERTEX WAS IN THE PLANE!\n");
		}
	}

	return ret;
}

Color sample(Texture *texture, float u, float v) {
	s32 x = (s32) (u * (float) (texture->w - 1));
	s32 y = (s32) (v * (float) (texture->h - 1));

	u32 pixel = ((u32 *) texture->data)[y * texture->w + x];

	u8 r = (u8) (pixel >> 0);
	u8 g = (u8) (pixel >> 8);
	u8 b = (u8) (pixel >> 16);
	u8 a = (u8) (pixel >> 24);

	Color color = {r, g, b, a};

	return color;
}

Vertex4 vertex_shader(Vertex v, Transform t) {
	glm::mat4x4 world_matrix = transform_get_matrix(&t);
	glm::mat4x4 camera_matrix = camera_get_matrix(&camera);

	glm::vec4 pos = camera_matrix * world_matrix * glm::vec4(v.pos, 1.0f);
	glm::vec3 norm = glm::mat3(glm::transpose(glm::inverse(world_matrix))) * v.norm;
	glm::vec2 tex = v.tex;

	return {pos, norm ,tex};
}

Color fragment_shader(s32 x, s32 y, glm::vec3 pos, glm::vec3 norm, glm::vec2 tex) {
	Color c = sample(&texture, tex.x, tex.y);

	return c;
}

void draw_triangle(Vertex v0, Vertex v1, Vertex v2) {
	glm::vec3 p0 = v0.pos;
	glm::vec3 p1 = v1.pos;
	glm::vec3 p2 = v2.pos;

	glm::vec3 n0 = v0.norm;
	glm::vec3 n1 = v1.norm;
	glm::vec3 n2 = v2.norm;

	glm::vec2 t0 = v0.tex;
	glm::vec2 t1 = v1.tex;
	glm::vec2 t2 = v2.tex;

	if (p0.y > p1.y) {
		std::swap(p0, p1);
		std::swap(n0, n1);
		std::swap(t0, t1);
	}

	if (p0.y > p2.y) {
		std::swap(p0, p2);
		std::swap(n0, n2);
		std::swap(t0, t2);
	}

	if (p1.y > p2.y) {
		std::swap(p1, p2);
		std::swap(n1, n2);
		std::swap(t1, t2);
	}

	// from y0 to y1

	s32 long_height = (s32) (p2.y - p0.y) + 1;

	for (s32 y = p0.y; y <= p1.y; ++y) {
		s32 short_height = (s32) (p1.y - p0.y) + 1; // + 1 so no div by zero
		float a = (float) ((float) y - (s32)p0.y) / (float) short_height;
		float b = (float) ((float) y - (s32)p0.y) / (float) long_height;

		s32 x_a = (s32) ((float) p0.x + a * (float) (p1.x - p0.x));
		s32 x_b = (s32) ((float) p0.x + b * (float) (p2.x - p0.x));

		glm::vec3 n_a = n0 + a * (n1 - n0);
		glm::vec3 n_b = n0 + b * (n2 - n0);

		glm::vec2 t_a = t0 + a * (t1 - t0);
		glm::vec2 t_b = t0 + b * (t2 - t0);

		float z_a = p0.z + a * (p1.z - p0.z);
		float z_b = p0.z + b * (p2.z - p0.z);

		if (x_a > x_b) {
			std::swap(x_b, x_a);
			std::swap(n_a, n_b);
			std::swap(t_a, t_b);
			std::swap(z_a, z_b);
		}

		for (s32 x = x_a; x <= x_b; ++x) {
			s32 length = x_b - x_a + 1;
			float t = (float) (x - x_a) / (float) length;

			glm::vec3 norm0 = n_a + t * (n_b - n_a);
			glm::vec2 tex0 = t_a + t * (t_b - t_a);

			float z = z_a + t * (z_b - z_a);

			glm::vec3 pos((float) x / z / (viewport_distance), (float) y / z / (viewport_distance), (float) z);
			glm::vec3 norm(norm0 / z);
			glm::vec2 tex(tex0 / z);

			Color c = fragment_shader(x, y, pos, norm, tex);

			set_pixel_center(x, y, z, c);
		}
	}

	// from y1 to y2

	for (s32 y = p1.y; y <= p2.y; ++y) {
		s32 short_height = (s32) (p2.y - p1.y) + 1; // + 1 so no div by zero
		float a = (float) (y - (s32)p1.y) / (float) short_height;
		float b = (float) (y - (s32)p0.y) / (float) long_height;

		s32 x_a = (int) ((float) p1.x + a * (float) (p2.x - p1.x));
		s32 x_b = (int) ((float) p0.x + b * (float) (p2.x - p0.x));

		glm::vec3 n_a = n1 + a * (n2 - n1);
		glm::vec3 n_b = n0 + b * (n2 - n0);

		glm::vec2 t_a = t1 + a * (t2 - t1);
		glm::vec2 t_b = t0 + b * (t2 - t0);

		float z_a = p1.z + a * (p2.z - p1.z);
		float z_b = p0.z + b * (p2.z - p0.z);

		if (x_a > x_b) {
			std::swap(x_b, x_a);
			std::swap(n_a, n_b);
			std::swap(t_a, t_b);
			std::swap(z_a, z_b);

		}

		for (s32 x = x_a; x <= x_b; ++x) {
			s32 length = x_b - x_a + 1;
			float t = (float) (x - x_a) / (float) length;

			glm::vec3 norm0 = n_a + t * (n_b - n_a);
			glm::vec2 tex0 = t_a + t * (t_b - t_a);

			float z = z_a + t * (z_b - z_a);

			glm::vec3 pos((float) x / z / (viewport_distance), (float) y / z / (viewport_distance), (float) z);
			glm::vec3 norm(norm0 / z);
			glm::vec2 tex(tex0 / z);

			Color c = fragment_shader(x, y, pos, norm, tex);

			set_pixel_center(x, y, z, c);
		}
	}
}

void draw_vertexlist_in_camera_space(Vertex4List vertexlist) {
	u32 vertex_count = vertexlist.count;

	assert(vertex_count % 3 == 0);

	if (vertex_count == 0) return;

	glm::vec3 projection_vertices[vertex_count];
	Vertex canvas_vertices[vertex_count];

	for (s32 i = 0; i < vertex_count; ++i) {
		projection_vertices[i] = projection_matrix * vertex4list_get(&vertexlist, i).pos;
		canvas_vertices[i].pos = canvas_matrix * projection_vertices[i];

		canvas_vertices[i].pos.x /= canvas_vertices[i].pos.z;
		canvas_vertices[i].pos.y /= canvas_vertices[i].pos.z;
		canvas_vertices[i].pos.z = 1.0f / vertex4list_get(&vertexlist, i).pos.z;

		canvas_vertices[i].norm = vertex4list_get(&vertexlist, i).norm / vertex4list_get(&vertexlist, i).pos.z;
		canvas_vertices[i].tex = vertex4list_get(&vertexlist, i).tex / vertex4list_get(&vertexlist, i).pos.z;
	}

	for (s32 i = 0; i < vertexlist.count - 2; i += 3) {

		Vertex4 p0 = vertex4list_get(&vertexlist, i);
		Vertex4 p1 = vertex4list_get(&vertexlist, i + 1);
		Vertex4 p2 = vertex4list_get(&vertexlist, i + 2);

		Vertex v0 = canvas_vertices[i];
		Vertex v1 = canvas_vertices[i + 1];
		Vertex v2 = canvas_vertices[i + 2];

		if (mode == TRIANGLE) {
			draw_triangle(v0, v1, v2);
		} else if (mode == LINE) {
			draw_triangle_wireframe(v0, v1, v2);
		}
	}
}

void draw_indexed(Vertex *vertices, u32 vertex_count, u32 *indices, u32 index_count, Transform t) {
	assert(index_count % 3 == 0);

	// local space
	// world transform
	// world space
	// camera transform
	// camera space
	// clipping
	// perspective transform
	// viewport space
	// canvas transform
	// canvas space (discrete)

	// first we take the whole vertex list through the vertex shader
	// then we cull the whole model if it is outisde the view frustum
	// then we do perspective projection for every vertex going from 3D homogenous to 2D homogenous, keeping the z coordinate of camera space in the 2D w coordinate of viewport space
	// then we assemble each primitive
	// then we clip the primitive
	// then we draw the primitive

	Vertex4 camera_vertices[vertex_count];

	for (s32 i = 0; i < vertex_count; ++i) {
		camera_vertices[i] = vertex_shader(vertices[i], t);
	}

	// culling the whole model

	// compute average of all the mesh's vertices
	glm::vec4 avg(0.0f, 0.0f, 0.0f, 0.0f);
	for (s32 j = 0; j < vertex_count; ++j) {
		avg += camera_vertices[j].pos;
	}
	avg /= vertex_count;

	float radius = 0.0f;
	for (s32 j = 0; j < vertex_count; ++j) {
		float distance = glm::length(camera_vertices[j].pos - avg);

		if (distance > radius) {
			radius = distance;
		}
	}

	if (signed_distance(near_plane, avg)   < -radius ||
		signed_distance(left_plane, avg)   < -radius ||
		signed_distance(right_plane, avg)  < -radius ||
		signed_distance(top_plane, avg)    < -radius ||
		signed_distance(bottom_plane, avg) < -radius) { return; }



	// whole model passed
	// primitive assembly
	// clipping
	// draw primitive

	for (int i = 0; i < index_count - 2; i += 3) {

		Vertex4 p0 = camera_vertices[indices[i]];
		Vertex4 p1 = camera_vertices[indices[i + 1]];
		Vertex4 p2 = camera_vertices[indices[i + 2]];

		// backface culling
		glm::vec3 v_1 = p1.pos - p0.pos;
		v_1 = glm::normalize(v_1);

		glm::vec3 v_2 = p2.pos - p0.pos;
		v_2 = glm::normalize(v_2);

		glm::vec3 n = glm::cross(v_1, v_2);

		glm::vec3 avg_v = (p0.pos + p1.pos + p2.pos) / 3.0f;

		glm::vec3 avg_to_camera = glm::vec3(0.0f, 0.0f, 0.0f) - avg_v;

		float result = glm::dot(avg_to_camera, n);

		if (result <= 0.0f) {
			printf("Culled triangle with indices %d, %d, %d becuase its normal was (%f, %f, %f)\n", indices[i], indices[i + 1], indices[i + 2], n.x, n.y, n.z);
//			debug_break();
			return;
		}

		Vertex4List list0 = {};
		vertex4list_init(&list0);

		vertex4list_add(&list0, p0);
		vertex4list_add(&list0, p1);
		vertex4list_add(&list0, p2);

		Vertex4List list1 = clip_vertexlist_against_plane(near_plane, list0);
		Vertex4List list2 = clip_vertexlist_against_plane(left_plane, list1);
		Vertex4List list3 = clip_vertexlist_against_plane(right_plane, list2);
		Vertex4List list4 = clip_vertexlist_against_plane(top_plane, list3);
		Vertex4List list5 = clip_vertexlist_against_plane(bottom_plane, list4);


		draw_vertexlist_in_camera_space(list5);

		vertex4list_free(&list0);
		vertex4list_free(&list1);
		vertex4list_free(&list2);
		vertex4list_free(&list3);
		vertex4list_free(&list4);
		vertex4list_free(&list5);
	}
}

int init() {

	SDL_Init(SDL_INIT_VIDEO);

	SDL_SetRelativeMouseMode(rel_mouse);

	counter_freq = SDL_GetPerformanceFrequency();

	window = SDL_CreateWindow("Rasterizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, window_width, window_height, SDL_WINDOW_RESIZABLE);

	backbuffer = SDL_GetWindowSurface(window);

	triplebuffer = SDL_CreateRGBSurfaceWithFormat(0, canvas_width, canvas_height, backbuffer->format->BitsPerPixel, backbuffer->format->format);

	depthbuffer = (float *) malloc(canvas_width * canvas_height * sizeof(float));

	float projection_array[12] = {
			viewport_distance, 0.0f,              0.0f,
			0.0f,              viewport_distance, 0.0f,
			0.0f,              0.0f,              1.0f,
			0.0f,              0.0f,              0.0f
	};

	projection_matrix = glm::make_mat4x3(projection_array);

	float canvas_array[9] = {
			(canvas_width - 1) / viewport_width, 0.0f,                                  0.0f,
			0.0f,                                (canvas_height - 1) / viewport_height, 0.0f,
			0.0f,                                0.0f,                                  1.0f
	};
	canvas_matrix = glm::make_mat3x3(canvas_array);

	camera_init(&camera);

	mode = TRIANGLE;

	texture.data = stbi_load("data/texture2.png", &texture.w, &texture.h, &texture.c, 4);

	keyboard = {};

	return 0;
}

void update(float dt) {
	if (keyboard.right) {
		camera.translation += (glm::vec3) {1.0f, 0.0f, 0.0f} * 5.0f * dt;
	}

	if (keyboard.left) {
		camera.translation += (glm::vec3) {-1.0f, 0.0f, 0.0f} * 5.0f * dt;
	}

	if (keyboard.forward) {
		camera.translation += (glm::vec3) {0.0f, 0.0f, 1.0f} * 5.0f * dt;
	}

	if (keyboard.backward) {
		camera.translation += (glm::vec3) {0.0f, 0.0f, -1.0f} * 5.0f * dt;
	}
}

void render() {
    clear_depthbuffer();
    clear((Color) {0xcc, 0xcc, 0xcc, 0x0});

    Model cube1 = {};

    Transform t1 = {
            glm::vec3(0.0f, 0, 5.0f),
            glm::vec3(glm::radians(0.0f), glm::radians(0.0f),glm::radians(0.0f)),
            glm::vec3(1.0f, 1.0f, 1.0f)
    };

    model_init(&cube1, cube_mesh, t1);

    model_render(&cube1);
}

int main(int argc, char *argv[]) {

	init();

	Vertex vertices[] = {
			{{1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {1.0f, 0.0f}},
			{{-1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 0.0f}},
			{{-1.0f, -1.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f}},
			{{1.0f, -1.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {1.0f, 1.0f}},
	};

	u32 indices[] = {
			0, 2, 1,
			0, 3, 2,
	};

	mesh_init(&cube_mesh, vertices, ARRAYCOUNT(vertices), indices, ARRAYCOUNT(indices));

	float old_time = (float) SDL_GetPerformanceCounter() / (float) counter_freq;

	while (running) {
		// Run the message loop.
		SDL_Event e;
		while (SDL_PollEvent(&e) != 0) {
			if (e.type == SDL_QUIT) {
				running = false;
			} else if (e.type == SDL_WINDOWEVENT) {
				switch (e.window.event) {
					case SDL_WINDOWEVENT_RESIZED: {
						window_width = e.window.data1;
						window_height = e.window.data2;

						SDL_Surface *new_backbuffer = SDL_GetWindowSurface(window);

						if (new_backbuffer != backbuffer) {
							backbuffer = new_backbuffer;
						}

						break;
					}
				}
			} else if (e.type == SDL_MOUSEMOTION) {
				camera.pitch += e.motion.yrel * 0.001f;
				camera.yaw += e.motion.xrel * 0.001f;
			} else if (e.type == SDL_KEYDOWN) {
				switch (e.key.keysym.sym) {
					case SDLK_ESCAPE: {
						SDL_Event qe;
						qe.type = SDL_QUIT;
						qe.quit.timestamp = SDL_GetTicks();
						SDL_PushEvent(&qe);
						break;
					}

					case SDLK_F1: {
						if (rel_mouse) {
							rel_mouse = SDL_FALSE;
							SDL_SetRelativeMouseMode(rel_mouse);
						} else {
							rel_mouse = SDL_TRUE;
							SDL_SetRelativeMouseMode(rel_mouse);
						}
						break;
					}

					case SDLK_F2: {
						if (rel_mouse) {
							rel_mouse = SDL_FALSE;
							SDL_SetRelativeMouseMode(rel_mouse);
						}

						debug_break();
						break;
					}

					case SDLK_F3: {
						if (mode == TRIANGLE) {
							mode = LINE;
						} else if (mode == LINE) {
							mode = TRIANGLE;
						}
						break;
					}

					case SDLK_d: {
						keyboard.right = true;
						break;
					}

					case SDLK_a: {
						keyboard.left = true;
						break;
					}

					case SDLK_w: {
						keyboard.forward = true;
						break;
					}

					case SDLK_s: {
						keyboard.backward = true;
						break;
					}
				}
			} else if (e.type == SDL_KEYUP) {
				switch (e.key.keysym.sym) {
					case SDLK_d: {
						keyboard.right = false;
						break;
					}

					case SDLK_a: {
						keyboard.left = false;
						break;
					}

					case SDLK_w: {
						keyboard.forward = false;
						break;
					}

					case SDLK_s: {
						keyboard.backward = false;
						break;
					}
				}
			}
		}

		float new_time = (float) SDL_GetPerformanceCounter() / (float) counter_freq;
		float frame_time = new_time - old_time;

		if (frame_time > 0.25f) frame_time = 0.25f;

		old_time = new_time;
		accumulator += frame_time;

		while (accumulator >= dt) {
            update(dt);
            t += dt;
            accumulator -= dt;
		}

		render();

		update_backbuffer();
		SDL_UpdateWindowSurface(window);

		printf("frame took %f s - %f fps\n", frame_time, 1.0f / frame_time);
		++frame;
	}

	SDL_DestroyWindow(window);
	SDL_Quit();

	return 0;
}
