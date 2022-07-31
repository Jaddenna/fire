use wasm_bindgen::{prelude::*, JsCast};
use web_sys::{Element, Document, EventListener, HtmlElement, console};
use std::fmt;
use rand::{Rng, prelude::SliceRandom};
use rand::thread_rng;
extern crate console_error_panic_hook;
use std::panic;

macro_rules! log {
    ( $( $t:tt )* ) => {
        web_sys::console::log_1(&format!( $( $t )* ).into());
    }
}

#[wasm_bindgen]
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cell {
    Dead = 0,
    Alive = 1,
}

#[wasm_bindgen]
pub struct Universe {
    width: u32,
    height: u32,
    cells: Vec<Cell>,
}

#[wasm_bindgen]
impl Universe {

    pub fn width(&self) -> u32 {
        self.width
    }

    pub fn height(&self) -> u32 {
        self.height
    }

    pub fn cells(&self) -> *const Cell {
        self.cells.as_ptr()
    }

    pub fn tick(&mut self) {
        let mut next = self.cells.clone();

        for row in 0..self.height {
            for column in 0..self.width {
                let index = self.get_index(row, column);
                let cell = self.cells[index];
                let live_nb = self.live_neighbor_count(row, column);

                let next_cell = match (cell, live_nb) {
                    (Cell::Alive, x) if x < 2 => Cell::Dead,
                    (Cell::Alive, 2) | (Cell::Alive, 3) => Cell::Alive,
                    (Cell::Alive, x) if x > 3 => Cell::Dead,
                    (Cell::Dead, 3) => Cell::Alive,
                    (otherwise, _) => otherwise, 
                };

                next[index] = next_cell;
            }
        }

        self.cells = next;
    }

    fn set_random_alive() -> Cell {
        if js_sys::Math::random() > 0.5 {
            return Cell::Alive;
        }        
        Cell::Dead
    }
    
    pub fn new() -> Universe {
        let width = 64;
        let height = 64;

        let cells = (0..width * height)
            .map(|_| {
                Self::set_random_alive()
            })
            .collect();

        Universe {
            width,
            height,
            cells,
        }
    }

    pub fn render(&self) -> String {
        self.to_string()
    }

    fn get_index(&self, row: u32, column: u32) -> usize {
        (row * self.width + column) as usize
    }

    fn live_neighbor_count(&self, row: u32, column: u32) -> u8 {
        let mut count = 0;
        for delta_row in [self.height - 1, 0 , 1].iter().cloned() {
            for delta_col in [self.width - 1, 0, 1].iter().cloned() {
                // eigene Position überspringen
                if delta_col == 0 && delta_row == 0 {
                    continue;
                }

                let nb_row = (row + delta_row) % self.height;
                let nb_col = (column + delta_col) % self.width;

                let index = self.get_index(nb_row, nb_col);

                count += self.cells[index] as u8;
            }
        }

        count
    }
}

impl fmt::Display for Universe {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for line in self.cells.as_slice().chunks(self.width as usize) {
            for &cell in line {
                let symbol = if cell == Cell::Dead { '◻' } else { '◼' };
                write!(f, "{}", symbol)?;
            }
            write!(f, "\n")?;
        }

        Ok(())
    }
}


pub struct Vector2 {
    x: f32,
    y: f32
}
 
impl Vector2 {   
    pub fn new(x: f32, y: f32) -> Vector2 {
        Vector2 { x: x, y: y }
    }
    
	pub fn dot(&self, other: Vector2) -> f32 {
		return self.x * other.x + self.y * other.y;
	}
}

pub fn shuffle(tab: Vec<u32>) -> Vec<u32> {
    let mut rng = rand::thread_rng();
    let mut next = tab.clone();

    for e in (next.len() - 1)..0 {
        let index = (rng.gen::<f32>() * (e as f32 - 1.0)).round() as usize;
        let temp = next[index];

        next[e] = tab[index];
        next[index] = temp;
    }

    next
}

pub fn make_permutation() -> Vec<u32> {
    let mut p:Vec<u32> = (0..256).collect();
    p.shuffle(&mut thread_rng());

    let mut copy = p.clone();
    p.append(&mut copy);
    p
}

pub fn get_constant_vector(v: u32) -> Vector2 {
	//v is the value from the permutation table
	let h = v & 3;

    let c_vec = match h {
        0 => Vector2::new(1.0, 1.0),
        1 => Vector2::new(-1.0, 1.0),
        2 => Vector2::new(-1.0, -1.0),
        _ => Vector2::new(1.0, -1.0),
    };

	c_vec
}

pub fn fade(t: f32) -> f32 {
	((6.0 * t - 15.0) * t + 10.0) * t * t * t
}

pub fn lerp(t: f32, a1: f32, a2: f32) -> f32 {
	a1 + t*(a2-a1)
}

pub fn noise_2d(P: Vec<u32>, x:f32, y: f32) -> f32{
    let x_floor = x.floor();
    let y_floor = y.floor();
    // log!("x {:?}, y {:?}", x, y);

	let X = (x_floor as u32 & 255) as usize;
	let Y = (y_floor as u32 & 255) as usize;
  
	let xf = x - x_floor;
	let yf = y - y_floor;
    // log!("x {:?}, y {:?}, x {:?}, y {:?}", X, Y, xf, yf);

	let top_right = Vector2::new(xf as f32 - 1.0, yf as f32 - 1.0);
	let top_left = Vector2::new(xf as f32, yf as f32 - 1.0);
	let bottom_right = Vector2::new(xf as f32 - 1.0, yf as f32);
	let bottom_left = Vector2::new(xf as f32, yf as f32);
	
	//Select a value in the array for each of the 4 corners
	let value_top_right = P[P[X + 1] as usize + Y + 1];
	let value_top_left = P[P[X] as usize + Y + 1];
	let value_bottom_right = P[P[X + 1] as usize + Y];
	let value_bottom_left = P[P[X] as usize + Y];
	
	let dot_top_right = top_right.dot(get_constant_vector(value_top_right));
	let dot_top_left = top_left.dot(get_constant_vector(value_top_left));
	let dot_bottom_right = bottom_right.dot(get_constant_vector(value_bottom_right));
	let dot_bottom_left = bottom_left.dot(get_constant_vector(value_bottom_left));
	
	let u = fade(xf as f32);
	let v = fade(yf as f32);
	
	return lerp(u,
		lerp(v, dot_bottom_left, dot_top_left),
		lerp(v, dot_bottom_right, dot_top_right)
	);

}

pub struct FireSettings {
    width: usize,
    height: usize,
    light_threshold: f32,
}

pub struct CoolingMap {
    pixels: Vec<u32>
}

impl CoolingMap {
    pub fn new(settings: FireSettings) -> CoolingMap{
        let width = 400;
        let height = 200;

        let P = make_permutation();
        let lightThreshold = settings.light_threshold;
        let frequency = 0.05;
        let mut map: Vec<Vec<u32>> = vec![];
        let mut repeat = true;
        let mut y = 0.0;
        let mut ch = 0.0;

        while repeat {//&& y < 1000.0 {
            let row = (0..width)
                .into_iter()
                .map(|x| ((noise_2d(P.clone(), x as f32 * frequency, y * frequency) + 1.0) * 0.5 * lightThreshold).round() as u32)
                .collect::<Vec<u32>>();
        
            // log!("row {:?}", row);

            let all_identical = y > 0.0 && map[0]
                .iter()
                .enumerate()
                .map(|(index, value)| *value == row[index])
                .all(|x| x);
            
            // log!("y {:?}, all_identical {:?}", y, all_identical);

            if all_identical || y == 999.0 {
                repeat = false;
                ch = y;
            } else {
                map.push(row);
            }

            y += 1.0;
        }

        let pixels = map.concat();

        CoolingMap {
            pixels
        }
    }

    // fn get_index(&self, row: u32, column: u32) -> usize {
    //     (row as usize * self.width + column as usize) as usize
    // }
}

#[wasm_bindgen]
pub struct Fire {
    width: usize,
    height: usize,
    pixels: Vec<u32>,
    cooling_map: Vec<u32>,
}

#[wasm_bindgen]
impl Fire {

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn pixels(&self) -> *const u32 {
        // log!("pixels {:?}", self.pixels.len());
        self.pixels.as_ptr()
    }

    pub fn last_p(&self) -> u32 {
        // log!("pixels last {:?}", self.pixels.last().unwrap());
        *self.pixels.last().unwrap()
    }

    pub fn cooling_map(&self) -> *const u32 {
        self.cooling_map.as_ptr()
    }

    pub fn tick(&mut self) {
        let mut buffer = self.pixels.clone();

        // log!("pixels {:?}", self.pixels);

        for row in 5..self.height - 1 {
            for column in 1..self.width - 1 {
                let rowU = row as u32;
                let columnU = column as u32;
                let indext = self.get_index(rowU - 1, columnU);
                let indexr = self.get_index(rowU, columnU + 1);
                let indexb = self.get_index(rowU + 1, columnU);
                let indexl = self.get_index(rowU, columnU - 1);
                let index = self.get_index(rowU, columnU);
                let cell = buffer[index];
                
                let nt = buffer[indext];
                let nb = buffer[indexb];
                let nl = buffer[indexl];
                let nr = buffer[indexr];
                let mut avg = (nt + nb + nl + nr) / 4;
                let cooling = self.cooling_map[index];
                avg = if avg > cooling { avg - cooling as u32 } else { 0 };
                let index5 = self.get_index(rowU - 3, columnU);
                self.pixels[index5] = avg;
            }
        }

        let mut map_chunks = self.cooling_map
            .chunks(self.width as usize)
            .map(|f| f.to_vec())
            .collect::<Vec<Vec<u32>>>();

        map_chunks.rotate_right(1);
        map_chunks.rotate_right(1);
        map_chunks.rotate_right(1);
        map_chunks.rotate_right(1);
        map_chunks.rotate_right(1);
        
        self.cooling_map = map_chunks.concat();

        // log!("pixels {:?}", self.pixels);
    }

    pub fn render(&self) -> String {
        self.to_string()
    }
    
    pub fn new() -> Fire {
        console_error_panic_hook::set_once();
        let width = 400;
        let height = 200;

        let P = make_permutation();
        let lightThreshold = 12.0;
        let frequency = 0.05;
        let mut map: Vec<Vec<u32>> = vec![];
        let mut repeat = true;
        let mut y = 0.0;
        let mut ch = 0.0;

        while repeat {//&& y < 1000.0 {
            let row = (0..width)
                .into_iter()
                .map(|x| ((noise_2d(P.clone(), x as f32 * frequency, y * frequency) + 1.0) * 0.5 * lightThreshold).round() as u32)
                .collect::<Vec<u32>>();
        
            // log!("row {:?}", row);

            let all_identical = y > 0.0 && map[0]
                .iter()
                .enumerate()
                .map(|(index, value)| *value == row[index])
                .all(|x| x);
            
            // log!("y {:?}, all_identical {:?}", y, all_identical);

            if all_identical || y == 999.0 {
                repeat = false;
                ch = y;
            } else {
                map.push(row);
            }

            y += 1.0;
        }
        // log!("ch {:?}", ch);

        let pixels: Vec<u32> = (0..(width * height))
            .enumerate()
            .map(|(index, value)| if index > (width * height) - width * 5 { 255 } else { 255 })
            .collect::<Vec<u32>>();

        let cooling_map = map.concat();
        // log!("y {:?}", pixels);
        Fire {
            width,
            height,
            pixels,
            cooling_map
        }
    }

    fn get_index(&self, row: u32, column: u32) -> usize {
        (row as usize * self.width + column as usize) as usize
    }
}

impl fmt::Display for Fire {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for line in self.cooling_map.as_slice().chunks(self.width as usize) {
            // for row in line {
            for &cell in line {
                write!(f, "|{}|", cell)?;
            }
            // }
            write!(f, "\n")?;
        }

        Ok(())
    }
}